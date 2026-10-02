import glob
import os, shutil
import tempfile
import pysam
import time
import random
import array
import multiprocessing


_COMPLEMENT_TABLE = str.maketrans({
    "A": "T",
    "C": "G",
    "T": "A",
    "G": "C",
})

_ION_P1B_3P_ADAPTER = "ATCACCGACTGCCCATAGAGAGGCTGAGAC"
# Tail length matching Torrent Suite --trim-zm behavior
_DEFAULT_ZM_TAIL_FLOWS = 16
_INT16_MIN = -32768
_INT16_MAX = 32767
_FLOW_SIZE = 480

# How the real TS caller works
# max_flow = min(
#     num_flows,
#     base_to_flow[last_base] + 16
# );
# 
# ...
# bam.AddTag("ZM", flowgram);

def _clip_int16(value):
    if value > _INT16_MAX:
        return _INT16_MAX
    if value < _INT16_MIN:
        return _INT16_MIN
    return int(value)

def remede(args, reporter):
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    if args.logfile:
      logfile_path = os.path.join(args.outdir, args.logfile)
    else:
       current_time = time.strftime("%Y%m%d-%H%M%S")
       logfile_path = os.path.join(args.outdir, f"{args.command}-{current_time}.log")
    reporter.set_logfile(logfile_path)

    reporter.info("Starting remede process...", color='c')
    bam_basename = os.path.splitext(os.path.basename(args.bam))[0]
    new_bam_path = os.path.join(args.outdir, f"reflow.{bam_basename}.bam")

    temp_bam_path = tempfile.mktemp(prefix=".temp_reflow_", suffix=".bam", dir=args.outdir)


    input_bam = pysam.AlignmentFile(args.bam, "rb")    
    reporter.info("Calculating total read counts...", color='w')
    read_count = input_bam.count()
    input_bam.close()
    input_bam = pysam.AlignmentFile(args.bam, "rb")
    # Write the intermediate BAM uncompressed and let the parallel (threaded)
    # samtools sort do the BGZF compression: compressed pysam writes run at
    # ~13k reads/s versus ~105k reads/s for uncompressed, and samtools sort -@
    # then compresses the whole file in a fraction of the time.
    output_bam = pysam.AlignmentFile(temp_bam_path, "wbu", template=input_bam)
    
    read_groups = input_bam.header["RG"]
    lib_key = "TCAG"
    barcode = "CAGATCCATCGAT"
    flow_order = "TACGTACGTCTGAGCATCGATCGATGTACAGC"
    processed_reads = 0
    reporter.info(f"Processing total {read_count} reads...", color='w')

    # Reads whose flow signal must not be recomputed.
    excluded_qnames = {
        "LKW3F:06614:10382",
        "LKW3F:02204:00972",
        "LKW3F:08640:03232",
        "LKW3F:08481:04892",
        "AEGIX:06652:01539",
    }

    # Parallelise the flow computation. The seed is derived from the read name
    # and sequence, so results are independent of scheduling order and the
    # output is identical to the serial path. BAM I/O and AlignedSegment
    # mutation stay in this process.
    workers = None
    parallelism = int(getattr(args, "processes", 0) or 0)
    if parallelism <= 0:
        parallelism = os.cpu_count() or 1
    if parallelism > 1 and read_count >= 500:
        workers = multiprocessing.Pool(processes=parallelism)

    adapter = _ION_P1B_3P_ADAPTER
    # Reads are streamed in chunks: one pool task computes the flow for a whole
    # chunk, so per-read IPC overhead is amortised. Chunk results are consumed
    # in stream order, keeping the output identical to the serial path.
    CHUNK = 64

    def apply_results(sam_chunk, results):
        """Apply a chunk's (flow_signal, final_seq) results in order.

        Reads are carried as SAM strings because pysam reuses one
        AlignedSegment across loop iterations; the string is the snapshot.
        """
        for sam, result in zip(sam_chunk, results):
            read = pysam.AlignedSegment.fromstring(sam, input_bam.header)
            flow_signal, final_seq = result
            read.set_tag("ZM", flow_signal, replace=True)
            if final_seq != read.query_sequence:
                qual = read.query_qualities
                final_len = len(final_seq)
                if qual is not None:
                    if final_len == 0:
                        new_qual = qual[:0]
                    elif read.is_reverse:
                        new_qual = qual[-final_len:]
                    else:
                        new_qual = qual[:final_len]
                else:
                    new_qual = None
                read.query_sequence = final_seq
                if new_qual is not None:
                    read.query_qualities = new_qual
            if read.cigartuples is None:
                continue
            trimmed_new_read = trim_cigar_tuples(read)
            extended_new_read = extend_cigar_tuples(trimmed_new_read)

            # Patched 2026-08-21 by Joowon Jang
            # Sum only the CIGAR operations that consume the query
            QUERY_CONSUMING_OPS = {0, 1, 4, 7, 8}  # M, I, S, =, X
            cigar_length = sum(
                length
                for operation, length in extended_new_read.cigartuples
                if operation in QUERY_CONSUMING_OPS
            )
            if extended_new_read.query_length != cigar_length:
                continue
            output_bam.write(extended_new_read)

    chunk_reads = []
    chunk_payloads = []
    inflight = []  # (reads_chunk, async_result)

    def dispatch_chunk():
        if not chunk_reads:
            return
        inflight.append((list(chunk_reads), workers.apply_async(_flow_chunk, (list(chunk_payloads),))))
        chunk_reads.clear()
        chunk_payloads.clear()

    def drain_inflight():
        """Consume at most one completed chunk, in submission order."""
        if not inflight:
            return
        reads_chunk, async_result = inflight.pop(0)
        apply_results(reads_chunk, async_result.get())

    for read in input_bam:
        processed_reads += 1
        if processed_reads % 1000 == 0:
            progress = (processed_reads / read_count) * 100
            reporter.progress(f"Progress: {processed_reads}/{read_count} ({progress:.2f}%)")
        elif processed_reads == read_count:
            reporter.progress(f"Progress: {processed_reads}/{read_count} (100.00%)")

        if read.qname in excluded_qnames:
            continue

        try:
            orig_flow_signal = read.get_tag("ZM")
        except KeyError:
            continue
        if not orig_flow_signal:
            continue

        sequence = read.query_sequence
        chunk_reads.append(read.to_string())
        chunk_payloads.append((
            read.qname,
            sequence,
            tuple(orig_flow_signal),
            read.is_reverse,
            flow_order,
            lib_key,
            barcode,
            adapter,
            None,
        ))

        if workers is not None:
            if len(chunk_reads) >= CHUNK:
                dispatch_chunk()
            if len(inflight) >= max(parallelism, 2):
                drain_inflight()
        else:
            apply_results([chunk_reads[-1]], [_flow_worker(chunk_payloads[-1])])
            chunk_reads.clear()
            chunk_payloads.clear()

    if workers is not None:
        dispatch_chunk()
        while inflight:
            drain_inflight()
        workers.close()
        workers.join()
    output_bam.close()

    shutil.copy(temp_bam_path, os.path.join(args.outdir, f"tmpreflow.{bam_basename}.bam"))

    # Step 5: Sort and index the new BAM file
    reporter.info("Sorting the output BAM file...", color='c')
    pysam.sort("-@", str(os.cpu_count()), "-o", new_bam_path, temp_bam_path)
    os.remove(temp_bam_path)

    reporter.info("Indexing the sorted BAM file...", color='c')
    pysam.index("-@", str(os.cpu_count()), new_bam_path)

    reporter.info("Cleaning up...", color='d')
    pattern = os.path.join(args.outdir, '.temp*')
    files_to_delete = glob.glob(pattern)
    for file_path in files_to_delete:
        if os.path.isfile(file_path):
            os.remove(file_path)

    reporter.info(f"Remede process completed. Output BAM: {new_bam_path}", color='g')

def _sequence_to_flow_counts(sequence, flow_order, max_flows):
    """
    Convert a sequence into Ion Torrent flow space.

    Returns:
        flow_counts:
            Number of bases incorporated at each flow.
            Example:
                sequence   = TCAG
                flow_order = TACG...

                -> [1, 0, 1, 0, 0, 1, 0, 1]

        seq_pos:
            Number of sequence bases consumed after processing max_flows.

    flow_order may be the full run flow order or a short repeat motif
    (e.g. TACG); modulo makes both work.
    """
    sequence = sequence.upper()
    flow_order = flow_order.upper()

    if not flow_order:
        raise ValueError("flow_order must not be empty")

    flow_order_len = len(flow_order)

    flow_counts = []
    seq_pos = 0
    flow_idx = 0
    seq_len = len(sequence)

    while seq_pos < seq_len and flow_idx < max_flows:
        flow_base = flow_order[flow_idx % flow_order_len]

        hp = 0
        while seq_pos < seq_len and sequence[seq_pos] == flow_base:
            hp += 1
            seq_pos += 1

        flow_counts.append(hp)
        flow_idx += 1

    return flow_counts, seq_pos

def _get_key_flow_info(
    lib_key,
    flow_order,
    orig_flow_signal,
    max_flows,
    default_one_mer_signal=256,
):
    """
    Compute the flow layout of the Ion Torrent library key.

    usable_key_flows = num_key_flows - 1

    The last key flow can form a homopolymer with the first nucleotide of
    the barcode/template, so the original ZM is not reused there.
    """
    key_flow_counts, consumed = _sequence_to_flow_counts(
        lib_key,
        flow_order,
        max_flows,
    )

    if consumed != len(lib_key):
        raise ValueError(
            f"lib_key cannot be represented within {max_flows} flows: "
            f"{lib_key!r}"
        )

    num_key_flows = len(key_flow_counts)
    usable_key_flows = max(0, num_key_flows - 1)

    # Estimate this read's 1-mer signal range from unambiguous key 1-mer flows
    one_mer_signals = []

    limit = min(
        usable_key_flows,
        len(orig_flow_signal),
    )

    for flow_idx in range(limit):
        if key_flow_counts[flow_idx] == 1:
            signal = int(orig_flow_signal[flow_idx])

            # Only normal 1-mer amplitudes feed the calibration
            if signal > 0:
                one_mer_signals.append(signal)

    if one_mer_signals:
        one_mer_min = min(one_mer_signals)
        one_mer_max = max(one_mer_signals)
    else:
        # ZM ~= normalized signal * 256, so 256 is used as the 1-mer fallback
        one_mer_min = default_one_mer_signal
        one_mer_max = default_one_mer_signal

    return (
        key_flow_counts,
        usable_key_flows,
        one_mer_min,
        one_mer_max,
    )

def _flow_worker(payload):
    """Compute the ZM signal and sequence for one read.

    The seed is derived from the read name and sequence, so the result is
    independent of scheduling order and identical to the serial path.
    """
    qname, sequence, orig_flow_signal, is_reverse, flow_order, lib_key, barcode, adapter, max_flows = payload
    return estimate_flow(
        sequence,
        orig_flow_signal,
        flow_order,
        lib_key,
        barcode,
        is_reverse,
        seed=f"{qname}\0{sequence}",
        adapter=adapter,
        max_flows=max_flows,
    )


def _flow_chunk(payloads):
    """Pool task: compute the flow for a whole chunk of reads."""
    return [_flow_worker(payload) for payload in payloads]


def estimate_flow(
    sequence,
    origFlowSignal,
    flowOrder,
    libKey,
    barcode,
    reverse,
    seed=None,
    *,
    adapter=_ION_P1B_3P_ADAPTER,
    max_flows=None,
    zm_tail_flows=_DEFAULT_ZM_TAIL_FLOWS,
    noise_min=-30,
    noise_max=30,
):
    """
    Reconstruct the Ion Torrent ZM flow signal.

    The template:

        libKey + barcode + read + 3' adapter

    is treated as one continuous DNA template.

    The original ZM is kept as-is up to usableKeyFlows; from the last key
    flow onward the signal is generated from the new template.
    """

    if not origFlowSignal:
        raise ValueError("origFlowSignal is empty")

    if not flowOrder:
        raise ValueError("flowOrder is empty")

    rng = random.Random(seed)

    # randint(a, b) == a + _randbelow(b - a + 1); _randbelow(n) draws
    # getrandbits(n.bit_length()) and rejects draws >= n (CPython random.py).
    # Inlining removes ~6 nested Python calls per draw while producing a
    # bit-identical stream, so output is unchanged.
    getrandbits = rng.getrandbits

    sequence = sequence.upper()
    libKey = libKey.upper()
    barcode = barcode.upper()
    adapter = adapter.upper()
    flowOrder = flowOrder.upper()

    # If the BAM is aligned to the reverse strand, restore the actual
    # sequencing orientation before computing the flows
    if reverse:
        oriented_sequence = reverse_complement(sequence)
    else:
        oriented_sequence = sequence

    #
    # Actual sequencing template
    #
    # key | barcode | insert | adapter
    #
    template = (
        libKey
        + barcode
        + oriented_sequence
        + adapter
    )

    prefix_len = len(libKey) + len(barcode)

    # Immediately before the adapter starts
    read_end = prefix_len + len(oriented_sequence)

    #
    # Total number of flows in the run
    #
    # When the Ion BAM's @RG FO holds the full flow string, len(flowOrder)
    # is the real num_flows.
    #
    # When only a short TACG motif is passed, the existing ZM length is used.
    #
    if max_flows is None:
        max_flows = max(
            len(flowOrder),
            len(origFlowSignal),
        )

    if max_flows <= 0:
        raise ValueError("max_flows must be positive")

    (
        key_flow_counts,
        usable_key_flows,
        one_mer_min,
        one_mer_max,
    ) = _get_key_flow_info(
        libKey,
        flowOrder,
        origFlowSignal,
        max_flows,
    )

    new_flow_signal = array.array("h")
    append_signal = new_flow_signal.append

    #
    # ------------------------------------------------------------
    # 1. Flows that involve only the key are preserved from the original ZM
    # ------------------------------------------------------------
    #
    preserve_flows = min(
        usable_key_flows,
        len(origFlowSignal),
        max_flows,
    )

    if preserve_flows:
        new_flow_signal.extend(
            origFlowSignal[:preserve_flows]
        )

    #
    # How many template bases have been consumed so far?
    #
    # Note:
    # usableKeyFlows excludes the last key flow, so for the usual TCAG
    # key this position sits just before the final G.
    #
    seq_pos = sum(
        key_flow_counts[:preserve_flows]
    )

    flow_order_len = len(flowOrder)
    template_len = len(template)

    #
    # Once the flow consuming the read's last base is known, extend the ZM by
    # +16 flows the way Torrent Suite does.
    #
    target_flow_count = None

    #
    # Special case: libKey/barcode/read already ended earlier
    #
    if seq_pos >= read_end:
        if preserve_flows:
            last_read_flow = preserve_flows - 1
            target_flow_count = min(
                max_flows,
                last_read_flow + zm_tail_flows,
            )
        else:
            target_flow_count = min(
                max_flows,
                zm_tail_flows,
            )

    #
    # ------------------------------------------------------------
    # 2. Synthetic flows computed from the last key flow onward
    # ------------------------------------------------------------
    #
    for flow_idx in range(preserve_flows, max_flows):

        flow_base = flowOrder[
            flow_idx % flow_order_len
        ]

        hp = 0

        #
        # The key -> barcode -> insert -> adapter boundaries are not
        # distinguished.
        #
        # Example:
        #
        #     ...G | GG...
        #
        # is treated as a 3-mer within the same G flow.
        #
        while (
            seq_pos < template_len
            and template[seq_pos] == flow_base
        ):
            hp += 1
            seq_pos += 1

        if hp:
            #
            # Use the original key's 1-mer signal range to generate the signal
            # of each incorporated nucleotide.
            #
            signal = 0

            # Inline randint(one_mer_min, one_mer_max): span/k are constant
            # here, so compute them once for the whole run of this base.
            # span >= 1 always (one_mer_min <= one_mer_max).
            span = one_mer_max - one_mer_min + 1
            k = span.bit_length()

            for _ in range(hp):
                draw = getrandbits(k)
                while draw >= span:
                    draw = getrandbits(k)
                signal += one_mer_min + draw

            append_signal(
                _clip_int16(signal)
            )

        else:
            #
            # 0-mer : randint(noise_min, noise_max), span = noise_max-noise_min+1
            #
            nspan = noise_max - noise_min + 1
            nk = nspan.bit_length()

            draw = getrandbits(nk)
            while draw >= nspan:
                draw = getrandbits(nk)

            append_signal(
                noise_min + draw
            )

        #
        # Did this flow consume the read's last base?
        #
        # Even when the adapter is consumed in the same homopolymer, the flow
        # holding the read's last base is the current flow.
        #
        if (
            target_flow_count is None
            and seq_pos >= read_end
        ):
            target_flow_count = min(
                max_flows,
                flow_idx + zm_tail_flows,
            )

        #
        # The required ZM tail has been generated
        #
        if (
            target_flow_count is not None
            and len(new_flow_signal) >= target_flow_count
        ):
            break

    #
    # ------------------------------------------------------------
    # 3. The actual BAM query_sequence does not include the adapter.
    # ------------------------------------------------------------
    #
    consumed_read_end = min(
        seq_pos,
        read_end,
    )

    final_oriented_sequence = template[
        prefix_len:consumed_read_end
    ]

    if reverse:
        final_seq = reverse_complement(
            final_oriented_sequence
        )
    else:
        final_seq = final_oriented_sequence

    return new_flow_signal, final_seq

def reverse_complement(seq):
    return seq.upper().translate(_COMPLEMENT_TABLE)[::-1]

def trim_cigar_tuples(read):
    query_length = read.query_length
    new_cigartuples = []
    cumulative_length = 0

    for operation, length in read.cigartuples:
        if operation in [0, 1, 4, 7, 8]:  # consider only M, I, S, =, X
            if cumulative_length + length <= query_length:
                new_cigartuples.append((operation, length))
                cumulative_length += length
            else:
                # Trim the overflowing part when there is one
                new_length = query_length - cumulative_length
                if new_length > 0:
                    new_cigartuples.append((operation, new_length))
                break
        else:
            new_cigartuples.append((operation, length))

    # Replace the CIGAR tuples
    read.cigartuples = new_cigartuples
    return read

def extend_cigar_tuples(read):
    query_length = read.query_length
    cigar_length = sum([length for (operation, length) in read.cigartuples if operation in [0, 1, 4, 7, 8]])

    if cigar_length < query_length:
        # Append an M (match) operation for the missing length
        missing_length = query_length - cigar_length
        read.cigartuples.append((0, missing_length))  # 0 means the M operation

    return read
