import glob
import os
import select
import shutil
import subprocess
import tempfile
import time
import sys
import re
import pysam
from reporting import safe_report

PROJECT_DIR = os.path.dirname(os.path.abspath(__file__))
BAMSURGEON_BIN = os.path.join(PROJECT_DIR, "deps", "bamsurgeon", "bin")
PICARD_JAR = os.path.join(PROJECT_DIR, "deps", "picard-tools-1.131", "picard.jar")

def check_java():
    p = subprocess.Popen(['java', '-version'], stderr=subprocess.PIPE)
    for line in p.stderr:
        line = line.decode()

        if line.startswith('java version') or line.startswith('openjdk version'):
            return True
    return False

def check_bwa():
    p = subprocess.Popen(['bwa'], stderr=subprocess.PIPE)
    for line in p.stderr:
        line = line.decode()

        if line.startswith('Version:'):
            major, minor, sub = line.strip().split()[1].split('.')
            sub = sub.split('-')[0]
            digit_pattern = re.compile(r'\D')
            sub = list(filter(None, digit_pattern.split(sub)))[0]
            if int(major) >= 0 and int(minor) >= 7 and int(sub) >= 12:
                return True
    return False

def check_samtools():
    p = subprocess.Popen(['samtools'], stderr=subprocess.PIPE)
    for line in p.stderr:
        line = line.decode()

        if line.startswith('Version:'):
            major, minor = line.strip().split()[1].split('.')[:2]
            minor = minor.split('-')[0]
            if int(major) >= 1 and int(minor) >= 2:
                return True
    return False

def check_wgsim():
    p = subprocess.Popen(['wgsim'], stderr=subprocess.PIPE)
    for line in p.stderr:
        line = line.decode()

        if line.startswith('Version:'):
            major, minor = line.strip().split()[1].split('.')[:2]
            minor = minor.split('-')[0]
            if int(major) >= 0 and int(minor) >= 2:
                return True
    return False

def check_velvet():
    p = subprocess.Popen(['velvetg'], stdout=subprocess.PIPE)
    for line in p.stdout: 
        line = line.decode()

        if line.startswith('Version'):
            major, minor = line.strip().split()[1].split('.')[:2]
            minor = minor.split('-')[0]
            if int(major) >= 1 and int(minor) >= 2:
                return True
    return False

def check_exonerate():
    p = subprocess.Popen(['exonerate'], stdout=subprocess.PIPE)
    for line in p.stdout:
        line = line.decode()

        if line.startswith('exonerate from exonerate'):
            major, minor = line.strip().split()[-1].split('.')[:2]
            minor = minor.split('-')[0]
            if int(major) >= 2 and int(minor) >= 2:
                return True
    return False

def check_python():
    return sys.version_info >= (3, 6)

def check_dependencies(reporter):
    if not check_python():
        safe_report(reporter, "error", 'Dependency problem: python >= 3.6 is required', color='r')
        return False
    if not check_bwa(): 
        safe_report(reporter, "error", 'Dependency problem: bwa >= 0.7.12 not found', color='r')
        return False
    if not check_samtools():
        safe_report(reporter, "error", 'Dependency problem: samtools >= 1.2 not found', color='r')
        return False
    if not check_wgsim():
        safe_report(reporter, "error", 'Dependency problem: wgsim not found (required for addsv)', color='r')
        return False
    if not check_velvet():
        safe_report(reporter, "error", 'Dependency problem: velvet >= 1.2 not found (required for addsv)', color='r')
        return False
    if not check_exonerate():
        safe_report(reporter, "error", 'Dependency problem: exonerate >= 2.2 not found (required for addsv)', color='r')
        return False
    if not check_java():
        safe_report(reporter, "error", 'Dependency problem: java not found', color='r')
        return False
    return True

def run_bamsurgeon(args, reporter, output_bam_path):
    bam_basename = os.path.splitext(os.path.basename(args.bam))[0]
    snv_bam_path = os.path.join(args.outdir, f"snv.{args.aligner}.{bam_basename}.bam")
    indel_bam_path = os.path.join(args.outdir, f"indel.{args.aligner}.{bam_basename}.bam")
    tmp_dir = os.path.join(args.outdir, "temp")
    tmp_snv_bam_path = tempfile.mktemp(prefix=".temp_", suffix=".bam", dir=args.outdir)    
    tmp_indel_bam_path = tempfile.mktemp(prefix=".temp_",suffix=".bam", dir=args.outdir)
    mindepth = int(getattr(args, "mindepth", 30))
    min_mutated_reads = int(getattr(args, "min_mutated_reads", 3))
    if mindepth < 2:
        raise ValueError("BAMSurgeon mindepth must be at least 2")
    if min_mutated_reads < 1:
        raise ValueError("BAMSurgeon min-mutated-reads must be at least 1")

    try:
        if args.snv:
            reporter.info("Fabricating SNVs...", color='c')
            reporter.info(f"Transfer the SNVs in {args.snv} to {args.bam}.", color='d')
            snv_cmd = []
            if args.aligner == "mem":
                snv_cmd = [
                    "python", "-O", os.path.join(BAMSURGEON_BIN, "addsnv.py"),
                    "-p", str(min(os.cpu_count(), 16)),
                    "-v", args.snv,
                    "-f", args.bam,
                    "-r", args.reference,
                    "-o", tmp_snv_bam_path,
                    "--picardjar", PICARD_JAR,
                    "--mindepth", str(mindepth),
                    "--aligner", args.aligner,
                    "--requirepaired",
                    "--tmpdir", tmp_dir,
                    "--vcf", os.path.join(args.outdir, "snv.vcf")
                ]
            elif args.aligner == "tmap":
                snv_cmd = [
                    "python", "-O", os.path.join(BAMSURGEON_BIN, "addsnv.py"),
                    "-p", str(min(os.cpu_count(),16)),
                    "-v", args.snv,
                    "-f", args.bam,
                    "-r", args.reference,
                    "-o", tmp_snv_bam_path,
                    "--picardjar", PICARD_JAR,
                    "--mindepth", str(mindepth),
                    "--aligner", args.aligner,
                    "--force",
                    "--insane",
                    "--single",
                    "--tagreads",
                    "--tmpdir", tmp_dir,
                    "--vcf", os.path.join(args.outdir, "snv.vcf")
                ]
            else:
                safe_report(reporter, "error", f"Unsupported aligner: {args.aligner}. Exiting.")
                return False
            if not getattr(args, "require_paired", True) and "--requirepaired" in snv_cmd:
                snv_cmd.remove("--requirepaired")
            if getattr(args, "force", False) and "--force" not in snv_cmd:
                snv_cmd.append("--force")
            if hasattr(args, "min_mutated_reads"):
                depth_index = snv_cmd.index("--mindepth") + 2
                snv_cmd[depth_index:depth_index] = [
                    "--minmutreads", str(min_mutated_reads)
                ]
            if getattr(args, "tag_reads", False) and "--tagreads" not in snv_cmd:
                snv_cmd.append("--tagreads")
            reporter.command(snv_cmd)
            pending = ["|", "/", "-", "\\"]
            proc = subprocess.Popen(
                snv_cmd, cwd=args.outdir, stdout=subprocess.PIPE,
                stderr=subprocess.PIPE, text=True
            )

            while True:
                ready_fds, _, _ = select.select([proc.stdout, proc.stderr], [], [], 0.1)  # Faster spinner update
                
                if proc.stdout in ready_fds:
                    output = proc.stdout.readline()
                    if output:
                        reporter.stdout(output.strip())
                
                if proc.stderr in ready_fds:
                    error = proc.stderr.readline()
                    if error:
                        reporter.stderr(error.strip())
                
                # Check if process has terminated
                if proc.poll() is not None:
                    break
                
                reporter.status(f"Mutagenesis... {pending[int(time.time() * 10) % 4]}")

            # Ensure all remaining output is read
            for output in proc.stdout:
                reporter.stdout(output.strip())
            for error in proc.stderr:
                reporter.stderr(error.strip())

        if not os.path.exists(tmp_snv_bam_path):
            safe_report(reporter, "error", "Fabricating SNV BAM was failed. Exiting.")
            return False
        
        reporter.info("Sorting the output BAM file...", color='c')
        pysam.sort("-@", str(os.cpu_count()), "-o", snv_bam_path, tmp_snv_bam_path)
        os.remove(tmp_snv_bam_path)
        reporter.info("Indexing the sorted BAM file...", color='c')
        pysam.index("-@", str(os.cpu_count()), snv_bam_path)    
        reporter.info(f"SNV BAM: {snv_bam_path}", color='g')

        if args.indel:
            reporter.info("Fabricating INDELs...", color='c')
            reporter.info(f"Transfer the INDELs in {args.indel} to {snv_bam_path}.", color='d')
            indel_cmd = []
            if args.aligner == "mem":
                indel_cmd = [
                    "python", "-O", os.path.join(BAMSURGEON_BIN, "addindel.py"),
                    "-p", str(min(os.cpu_count(),16)),
                    "-v", args.indel,
                    "-f", snv_bam_path,
                    "-r", args.reference,
                    "-o", tmp_indel_bam_path,
                    "--picardjar", PICARD_JAR,
                    "--mindepth", str(mindepth),
                    "--aligner", args.aligner,
                    "--requirepaired",
                    "--tmpdir", tmp_dir,
                    "--vcf", os.path.join(args.outdir, "indel.vcf")
                ]
            elif args.aligner == "tmap":        
                indel_cmd = [
                    "python", "-O", os.path.join(BAMSURGEON_BIN, "addindel.py"),
                    "-p", str(min(os.cpu_count(),16)),
                    "-v", args.indel,
                    "-f", snv_bam_path,
                    "-r", args.reference,
                    "-o", tmp_indel_bam_path,
                    "--picardjar", PICARD_JAR,
                    "--mindepth", str(mindepth),
                    "--aligner", args.aligner,
                    "--force",
                    "--insane",
                    "--single",
                    "--tagreads",
                    "--tmpdir", tmp_dir,
                    "--vcf", os.path.join(args.outdir, "indel.vcf")
                ]            
            else:
                safe_report(reporter, "error", f"Unsupported aligner: {args.aligner}. Exiting.")
                return False
            if not getattr(args, "require_paired", True) and "--requirepaired" in indel_cmd:
                indel_cmd.remove("--requirepaired")
            if getattr(args, "force", False) and "--force" not in indel_cmd:
                indel_cmd.append("--force")
            if hasattr(args, "min_mutated_reads"):
                depth_index = indel_cmd.index("--mindepth") + 2
                indel_cmd[depth_index:depth_index] = [
                    "--minmutreads", str(min_mutated_reads)
                ]
            if getattr(args, "tag_reads", False) and "--tagreads" not in indel_cmd:
                indel_cmd.append("--tagreads")
            reporter.command(indel_cmd)
            pending = ["|", "/", "-", "\\"]
            proc = subprocess.Popen(
                indel_cmd, cwd=args.outdir, stdout=subprocess.PIPE,
                stderr=subprocess.PIPE, text=True
            )

            while True:
                ready_fds, _, _ = select.select([proc.stdout, proc.stderr], [], [], 0.1)  # Faster spinner update
                
                if proc.stdout in ready_fds:
                    output = proc.stdout.readline()
                    if output:
                        reporter.stdout(output.strip())
                
                if proc.stderr in ready_fds:
                    error = proc.stderr.readline()
                    if error:
                        reporter.stderr(error.strip())
                
                # Check if process has terminated
                if proc.poll() is not None:
                    break
                
                reporter.status(f"Mutagenesis... {pending[int(time.time() * 10) % 4]}")

            # Ensure all remaining output is read
            for output in proc.stdout:
                reporter.stdout(output.strip())
            for error in proc.stderr:
                reporter.stderr(error.strip())

            #while True:
            #    output = proc.stderr.readline()
            #    if output == '' and proc.poll() is not None:
            #        break
            #    if output:
            #        reporter.stdout(output.strip())
        if not os.path.exists(tmp_indel_bam_path):
            if getattr(args, "allow_partial", False) and os.path.exists(snv_bam_path):
                safe_report(
                    reporter,
                    "warning",
                    "No INDEL mutation succeeded; preserving the SNV-mutated BAM so validation can mark missed variants.",
                )
                shutil.copy(snv_bam_path, output_bam_path)
                pysam.index("-@", str(os.cpu_count()), output_bam_path)
                reporter.info(f"Partially completed: {output_bam_path}", color='y')
                return True
            safe_report(reporter, "error", "Fabricating INDEL BAM was failed. Exiting.")
            return False
        
        reporter.info("Sorting the INDEL BAM file...", color='c')
        pysam.sort("-@", str(os.cpu_count()), "-o", indel_bam_path, tmp_indel_bam_path)
        os.remove(tmp_indel_bam_path)

        reporter.info("Indexing the sorted INDEL BAM file...", color='c')
        pysam.index("-@", str(os.cpu_count()), indel_bam_path)
        reporter.info(f"Indel BAM: {indel_bam_path}", color='g')
        shutil.copy(indel_bam_path, output_bam_path)
        pysam.index("-@", str(os.cpu_count()), output_bam_path)
        reporter.info(f"Completed: {output_bam_path}", color='g')
        return True
    finally:
        if os.path.exists(tmp_snv_bam_path):
            os.remove(tmp_snv_bam_path)
        if os.path.exists(tmp_indel_bam_path):
            os.remove(tmp_indel_bam_path)
        if os.path.exists(tmp_dir):
            shutil.rmtree(tmp_dir)

def oracode(args, reporter):
    try:
        if not os.path.exists(args.outdir):
            os.makedirs(args.outdir)

        if args.logfile:
            logfile_path = os.path.join(args.outdir, args.logfile)
        else:
            current_time = time.strftime("%Y%m%d-%H%M%S")
            logfile_path = os.path.join(args.outdir, f"{args.command}-{current_time}.log")
        reporter.set_logfile(logfile_path)

        reporter.info("Checking dependencies...", color='c')
        if not check_dependencies(reporter):
            safe_report(reporter, "error", "Dependencies are not satisfied. Exiting.")
            raise Exception("Dependencies are not satisfied. Exiting.")   
        
        output_bam_path = os.path.join(args.outdir, f"bs.{args.aligner}.{os.path.splitext(os.path.basename(args.bam))[0]}.bam")
        if not run_bamsurgeon(args, reporter, output_bam_path):
            safe_report(reporter, "error", "BAM files are not generated. Exiting.")
            raise Exception("BAM files are not generated. Exiting.")
    except Exception as e:
        safe_report(reporter, "error", f"{str(e)}")
        raise
    finally:
        # Clean up
        safe_report(reporter, "info", "Cleaning up...", color='d')
        pattern = os.path.join(args.outdir, '*.temp*')
        files_to_delete = glob.glob(pattern)
        for file_path in files_to_delete:
            try:
                if os.path.isfile(file_path):
                    os.remove(file_path)
                elif os.path.isdir(file_path):
                    shutil.rmtree(file_path)
            except OSError as exc:
                safe_report(reporter, "warning", f"Cleanup could not remove {file_path}: {exc}")

        pattern = os.path.join(args.outdir, 'tmp.*')
        files_to_delete = glob.glob(pattern)
        for file_path in files_to_delete:
            try:
                if os.path.isfile(file_path):
                    os.remove(file_path)
                elif os.path.isdir(file_path):
                    shutil.rmtree(file_path)
            except OSError as exc:
                safe_report(reporter, "warning", f"Cleanup could not remove {file_path}: {exc}")

        pattern = os.path.join(args.outdir, 'addsnv*')
        files_to_delete = glob.glob(pattern)
        for file_path in files_to_delete:
            try:
                if os.path.isfile(file_path):
                    os.remove(file_path)
                elif os.path.isdir(file_path):
                    shutil.rmtree(file_path)
            except OSError as exc:
                safe_report(reporter, "warning", f"Cleanup could not remove {file_path}: {exc}")

        pattern = os.path.join(args.outdir, 'addindel*')
        files_to_delete = glob.glob(pattern)
        for file_path in files_to_delete:
            try:
                if os.path.isfile(file_path):
                    os.remove(file_path)
                elif os.path.isdir(file_path):
                    shutil.rmtree(file_path)
            except OSError as exc:
                safe_report(reporter, "warning", f"Cleanup could not remove {file_path}: {exc}")
