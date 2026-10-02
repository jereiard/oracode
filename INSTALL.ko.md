# 설치

> English: [INSTALL.md](INSTALL.md)

Oracode는 이 저장소에서 빌드한 Docker 이미지로 실행됩니다.

## 요구사항

- Docker가 설치된 Linux 호스트 (런처가 현재 디렉터리와 Docker 소켓을 마운트합니다).
- 이미지 및 작업 파일용 여유 디스크 약 20 GB.
- GPU는 필수가 아닙니다. 런처는 선택적 DeepVariant 경로를 위해 `--gpus all`을
  전달하며, GPU가 없으면 해당 경로만 건너뜁니다.
- **Ion Torrent 툴체인 이미지(필수).** Ion Torrent 변이 호출은 `tmap`, `tvc`,
  `tvcutils`, `tvcassembly`를 사용하며, 이는 재배포가 불가합니다. Torrent Suite
  소스(https://github.com/iontorrent/TS)의 지침에 따라 직접 이미지를 빌드해야
  합니다(해당 저장소의 `buildTools/BUILD.txt` 참고, Analysis 모듈은
  `MODULES=Analysis ./buildTools/build.sh`로 빌드). 빌드 시
  `--build-arg TS_IMAGE=<tsuite-이미지>`로 지정하십시오(기본값
  `jereiard/tsuite:5.18.1`). 그 외 구성요소는 모두 공개 소스에서 빌드됩니다.

## 이미지 빌드

이 저장소에서 이미지를 직접 빌드합니다. Dockerfile에는 필요한 리소스가 포함되어 있습니다. base
(Ubuntu 20.04 + Miniforge + samtools/bcftools + BAMSurgeon + Picard)와 최종
이미지 모두 공개 소스에서 이 저장소 안에서 빌드됩니다. 외부 입력은 Ion Torrent
툴체인(`TS_IMAGE`)뿐입니다. base 패키지, BAMSurgeon, Picard, IGV를 내려받기 위해
네트워크가 필요하며, `dist.sh`가 프록시 환경변수를 전달합니다.

```bash
# 빌드 (HTTP(S)_PROXY / NO_PROXY가 설정되어 있으면 자동 전달)
./dist.sh

# 빌드된 태그로 실행
ORACODE_IMAGE=oracode ./oracode --help
```

빌드 인자:

| 인자 | 기본값 | 용도 |
|---|---|---|
| `TS_IMAGE` | `jereiard/tsuite:5.18.1` | `tmap`, `tvc`, `tvcutils`, `tvcassembly`, TVC 파라미터 세트 제공. Ion Torrent 변이 호출에 필수. |
| `BAMSURGEON_REPO` | `https://www.github.com/jereiard/bamsurgeon` | BAMSurgeon 소스. |
| `BAMSURGEON_COMMIT` | `50b14964…` | 고정된 BAMSurgeon 커밋. |
| `PICARD_VERSION` | `1.131` | Picard tools 릴리스. |
| `IGV_VERSION` / `IGV_SHA256` | `2.19.8` / 고정값 | IGV 다운로드 및 체크섬. |

빌드 검증:

```bash
docker run --rm oracode --help
docker run --rm oracode remede --help
```

## 실행

`./oracode`는 `docker run`을 감싸며 모든 인자를 그대로 전달합니다.

```bash
# 입력 검증 및 표적/BAMSurgeon 파일 생성만 수행
./oracode --ui plain proficiency materials/variants.csv --prepare-only

# 4개 프로파일 전체 실행 (hg19/hg38 × Illumina/IonTorrent)
./oracode proficiency materials/variants.csv

# 데이터 루트를 명시적으로 지정
./oracode proficiency materials/variants.csv --data-root /path/to/data-root
```

이미지 내부 셸 열기:

```bash
./console
```

## 데이터 배치

Data root에는 `materials/`, `sources/`, `references/`가 있어야 합니다.
정확한 디렉터리 구조(변이 CSV 컬럼, 참조 FASTA 인덱스, 플랫폼별 소스
BAM/BED 경로)는 [README.ko.md](README.ko.md)의 **입력** 절에 명시되어
있습니다.

## 폐쇄망(air-gapped) 호스트

IGV는 시작 시 Google OAuth 조회를 수행하며, 폐쇄망 호스트에서는 headless 배치
실행이 프록시 407 또는 타임아웃으로 중단됩니다. 이미지에는 offline OAuth 설정과
`igv-headless`가 포함되어 있고, `run_igv`가 이를 사용하도록 지정하므로 네트워크
호출이 발생하지 않습니다.

## 문제 해결

| 증상 | 원인 / 해결 |
|---|---|
| `tvc: command not found` | TS 도구가 복사되지 않았습니다. `TS_IMAGE` 빌드 인자를 확인하십시오. |
| IGV가 시작 시 멈춤 | OAuth 조회 문제입니다. 이미지에 `igv-headless`와 offline OAuth 설정이 있는지 확인하십시오. |
| 프록시 빌드 실패 | `HTTP_PROXY`/`HTTPS_PROXY`/`NO_PROXY`를 export한 뒤 `./dist.sh`를 사용하십시오. |
| 출력 BAM을 TVC가 읽지 못함 | flow signal이 재생성되지 않았습니다. 변이 호출 전에 `remede`가 실행되는지 확인하십시오. |
