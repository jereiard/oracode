# Oracode

> English: [README.md](README.md)

Oracode는 germline NGS 정도관리물질을 **in silico**로 제작하는 도구입니다.
공개 정렬 BAM에서 표적 영역을 추출하고, BAMSurgeon으로 정의된 변이를 도입한 뒤,
Ion Torrent flow signal을 재작성하여 편집된 read가 변이 호출기에 인식되도록 하고,
변이 호출(Ion Torrent는 TVC, Illumina는 DeepVariant)과 변이별 IGV 스냅샷/보고서를
생성합니다.

핵심 단계는 **Ion Torrent flow signal 재생성**입니다. query sequence만 수정하고
ZM/FO/KS 태그를 갱신하지 않은 BAM은 Torrent Variant Caller가 해석하지 못하므로
도입한 변이가 호출되지 않습니다. Oracode는 flow order와 read별 결정적 seed로
ZM 배열을 재계산합니다.

## 구성

| 경로 | 용도 |
|---|---|
| `reflow.py` | CLI 진입점 (`proficiency`, `extractROI`, `bespoke`, `oracode`, `remede`, `educeVariant`, `igvbatch`) |
| `proficiency.py` | 4개 프로파일 전체 workflow (hg19/hg38 × Illumina/IonTorrent) |
| `oracode.py` | 서열 수준 변이 도입 (BAMSurgeon) |
| `remede.py` | Ion Torrent flow signal(ZM) 재생성 |
| `bespoke.py` | 원본 BAM에서 관심 영역(ROI) 추출 |
| `educe_variant.py` | 변이 호출(TVC / DeepVariant) 및 VCF 병합 |
| `igvbatch.py`, `igv-headless.sh` | headless IGV 스냅샷 생성 |
| `extract_roi.py` | gene list 기반 표적 영역 정의 |
| `reporting.py`, `textual_ui.py`, `curses_utils.py`, `logfile_manager.py` | 콘솔 리포팅 및 UI |
| `Dockerfile` | 이미지 빌드 (base + 최종 이미지) |
| `conda-explicit.txt` | base용 conda 패키지 고정 목록 |
| `patches/` | BAMSurgeon 패치 |
| `oracode`, `console`, `dist.sh` | 호스트 런처 및 이미지 빌드 스크립트 |

## 설치

[INSTALL.ko.md](INSTALL.ko.md)를 참고하십시오. 이 저장소에서 이미지를 직접
빌드합니다. 외부 입력은 Ion Torrent 툴체인(`tmap`/`tvc`)뿐이며,
[Torrent Suite 소스](https://github.com/iontorrent/TS)에서 빌드한 이미지를
`TS_IMAGE`로 지정합니다.

## 빠른 시작

```bash
./oracode --ui plain proficiency materials/variants.csv --prepare-only
./oracode proficiency materials/variants.csv
```

## 입력

입력은 변이 CSV와 data root 두 가지입니다. data root는 `materials/`,
`sources/`, `references/`를 포함하는 디렉터리로 자동 인식되며,
`--data-root`로 명시할 수 있습니다.

### 변이 CSV

다음 컬럼을 가진 헤더 행이 필요합니다(예시: `materials/variants.csv`).

| 컬럼 | 예 | 의미 |
|---|---|---|
| `contig` | `chr17` | 염색체 |
| `pos_hg19` | `7123443` | hg19 기준 1-기반 위치 |
| `pos_hg38` | `7220124` | hg38 기준 1-기반 위치 |
| `ref` | `C` | 참조 대립유전자 |
| `alt` | `A` | 대체 대립유전자 |
| `gene` | `ACADVL` | 유전자 심볼 |
| `type` | `SNV` | `SNV`, `INS`, `DEL` |

다음 두 컬럼은 선택 사항입니다.

| 컬럼 | 예 | 의미 |
|---|---|---|
| `vaf` | `0.25` | 도입할 변이 대립유전자 비율. `0 < vaf <= 1`. 셀이 비어 있거나 컬럼 자체가 없으면 `--vaf`(기본값 `0.5`)를 사용합니다. |
| `id` | `CVSET1-042` | `variants.results.csv`에 기록되는 변이 식별자. 컬럼이 없거나 셀이 비면 좌위에서 `gene:contig:pos_hg19:ref>alt`로 생성되며, 같은 좌위를 다시 나열한 행에는 `#2`, `#3` … 이 붙습니다(예: `ACADVL:chr17:7123443:C>A#2`). 직접 지정시에는 식별자는 고유하여야 합니다. |

`vaf`는 변이별 BAMSurgeon 입력 파일에 기록되므로, 전역 단일 비율이 아니라
변이마다 지정한 비율로 도입됩니다. 범위를 벗어나거나 숫자가 아닌 값은 실행을
중단시킵니다.

### Data root 구조

```
materials/
  variants.csv                              # 기본 입력 (경로는 임의 가능)
  gencode.v49lift37.annotation.sorted.gtf.gz   # hg19 표적 (proficiency)
  gencode.v49.annotation.sorted.gtf.gz         # hg38 표적 (proficiency)
  <gene list>                               # extractROI 전용 (옵션)
  <refseq>.db                               # extractROI 전용 (옵션)
references/
  hg19.chrs.fa  (+ .fai, .bwt/.pac/.sa/.amb/.ann, .tmap.*)   # bwa + tmap 인덱스
  hg38.chrs.fa  (+ .fai, .bwt/.pac/.sa/.amb/.ann, .tmap.*)
sources/
  illumina/hg19/HG001/HG001.bam  (+ .bai)
  illumina/hg19/HG001/HG001.bed
  illumina/hg38/HG001/HG001.bam  (+ .bai)
  illumina/hg38/HG001/HG001.bed
  iontorrent/hg19/HG001/HG001.bam (+ .bai)
  iontorrent/hg19/HG001/HG001.bed
  iontorrent/hg38/HG001/HG001.bam (+ .bai)
  iontorrent/hg38/HG001/HG001.bed
```

`proficiency`가 읽는 것은 GENCODE GTF(`materials/`)와 소스 BAM/BED
(`sources/`)뿐입니다. gene list와 RefSeq DB는 **`proficiency`에서 사용되지
않으며**, 별도 서브커맨드 `extractROI`에 속합니다.

참고:

- 참조 FASTA는 **bwa 인덱스**(Illumina `--aligner mem`)와 **tmap 인덱스**
  (Ion Torrent `--aligner tmap`)를 모두 갖추어야 합니다.
- `sources/<platform>/<genome>/HG001/HG001.bed`는 해당 플랫폼의 panel/ROI
  정의입니다.
- Ion Torrent 소스 BAM은 `ZM`/`FO`/`KS` 태그를 지닌 Torrent Suite 재정렬
  BAM이어야 합니다. Illumina 소스 BAM은 표준 paired-end BAM입니다.
- `--profiles`로 지정한 프로파일의 파일만 있으면 됩니다. 예를 들어
  `--profiles iontorrent-hg38`은 해당 소스 디렉터리 하나만 필요합니다.

## 라이선스

Apache License 2.0. [LICENSE](LICENSE)를 참고하십시오.
