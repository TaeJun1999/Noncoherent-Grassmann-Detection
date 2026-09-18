# 프로젝트 규칙

## 환경

- **실행 위치**: Docker 컨테이너 (`9000babdde68`), 계정 `HTJ`, 사용자별 분리 컨테이너
- **시각**: `TZ=Asia/Seoul`. `date`와 로그 타임스탬프는 KST (UTC 아님)
- **conda**: Miniforge (`~/miniforge3`), 채널은 conda-forge만. `defaults` 채널 사용 금지
  - `base` (Python 3.14): 도구 전용(tmux, node, uv). **실험 패키지 설치 금지**
  - `torch` (Python 3.12, torch 2.14.0+cu130): 실험용. 모든 파이썬 실행은 이 환경
- **파이썬 실행 방식**: 활성 환경을 가정하지 말고 절대경로로 호출한다
  - 실행: `~/miniforge3/envs/torch/bin/python …`
  - 설치: `~/miniforge3/envs/torch/bin/python -m pip install …`
  - 사람이 터미널에서 쓸 때는 `conda activate torch`로 충분
- **작업 디렉토리**: `~/ng`. Claude Code는 반드시 여기서 실행 (홈에서 띄우면 프로젝트 설정과 사용자 설정 파일이 겹침)
- **데이터 경로**: `~/data` (공용 경로 없음, 홈 내부). `~/ng/data`는 여기로 가는 심볼릭 링크
- **도구 경로**: `~/.local/bin`에 claude, node, tmux, uv, graphify가 있어 어느 conda 환경에서든 보인다

## 설치 규칙

- **apt / yum 등 시스템 패키지 설치 금지.** sudo 권한이 있어도 사용하지 않는다
- 홈 디렉토리(`~`) 밖에 파일을 만들지 않는다 — 컨테이너 레이어라 재생성 시 소멸한다
- 실험 패키지는 `torch` 환경 안에서만 pip(또는 conda)로 설치. base와 시스템 파이썬에는 설치 금지
- PyTorch는 conda가 아니라 pip로: `pip install torch torchvision --index-url https://download.pytorch.org/whl/cu130`
- 새 CLI 도구가 필요하면: 1) `conda install -c conda-forge <도구>` (base에 들어감) 2) 파이썬 기반 CLI는 `uv tool install <패키지>` (`~/.local/bin`에 들어감)

## GPU

- **하드웨어**: NVIDIA RTX PRO 6000 Blackwell Server Edition × 6, 각 96GB
- **CUDA 13.0 / 드라이버 580.173.02 / sm_120 (capability 12.0)**
- **Compute Mode: Exclusive_Process** — GPU당 프로세스 1개만 가능

규칙:
- 학습·추론 실행 전 반드시 `nvidia-smi`로 빈 GPU 확인
- `CUDA_VISIBLE_DEVICES`를 항상 명시 (예: `CUDA_VISIBLE_DEVICES=0 ~/miniforge3/envs/torch/bin/python src/train.py`)
- 점유된 GPU 접근 시 에러가 나는 것이 정상 동작이다. 자동 재시도를 반복하지 말고 사용자에게 알릴 것
- 장시간 학습은 실행 전 사용자 승인을 받을 것
- PyTorch를 새로 설치하면 `torch.cuda.get_arch_list()`에 `sm_120`이 있는지 확인. 오래된 빌드는 커널이 맞지 않는다 (2026-09-16 검증 완료: torch 2.14.0+cu130)

## 알려진 제약: /dev/shm 64MB

Docker 기본값이라 PyTorch DataLoader의 공유 메모리가 부족하다.

```
RuntimeError: DataLoader worker (pid xxx) is killed by signal: Bus error
```

- `num_workers`는 **2 이하**로 유지. 기본은 `0`
- 이 에러가 나면 워커 수 문제이므로 다른 원인을 찾지 말 것
- 관리자가 `--shm-size`를 늘려주면 이 제약을 해제하고 이 항목을 수정한다

## 자원 사용

서버를 다른 사용자와 공유하므로:
- `num_workers` 상향 금지 (위 shm 제약과 별개로 CPU 부담)
- 대용량 전처리나 고부하 작업은 실행 전 사용자에게 확인
- TensorBoard / Jupyter 기본 포트(6006, 8888) 사용 금지. 별도 포트 지정
- 오래된 체크포인트는 주기적으로 정리. 디스크는 공유 자원이다

## 실험 관리

- **하이퍼파라미터는 `configs/`의 YAML로만 관리.** 코드에 하드코딩 금지
- 실행 로그는 파일로 남길 것: `… 2>&1 | tee logs/$(date +%Y%m%d_%H%M).log`
- 긴 학습은 tmux 안에서 실행
- 커밋 메시지에 관련 실험 설정을 명시하면 나중에 추적이 쉽다
- 실험이 끝나면 즉시 `docs/EXPERIMENTS.md`에 한 행 추가: 날짜, 설정 파일, seed, 커밋 해시, GPU, 핵심 지표, 로그 경로, 메모. 이 표는 실행 사실의 기록이고, 인용 가능한 숫자는 `docs/02_RESULTS_LOG.md`에 있다
- 실험 결과를 문서에 반영할 때는 `docs/02_RESULTS_LOG.md`에 항목을 추가한다. 로그 원문의 숫자만 옮기고, 해석과 결론은 쓰지 않는다 — 그건 저자가 웹에서 한다
- 웹에서 만든 계획 문서는 `docs/plans/`에 들어온다. 구현 요청에 계획 파일이 지정되면 그것을 먼저 읽는다

## Git

- `data/`, `checkpoints/`, `outputs/`, `wandb/`, `logs/`, `*.pt`, `*.ckpt`, `graphify-out/`, `.claude/settings.local.json`은 커밋 금지 (.gitignore에 등록됨)
- 의미 있는 변경 단위로 커밋. 디버깅 중 임시 수정은 모아서 정리 후 커밋

## 코드 스타일

- 새 의존성을 추가하기 전에 이미 설치된 것으로 해결되는지 확인
- 실험 코드라도 재현 가능해야 한다: seed 고정, 설정 저장, 버전 기록
- 단, ablation 분기·로깅·중간 결과 저장은 "불필요한 코드"가 아니다. 연구 코드의 필수 요소다

## 작업 방식

- 큰 변경 전에 계획을 먼저 제시할 것
- 파일을 수정하기 전에 현재 내용을 읽고 확인할 것
- 에러가 나면 추측으로 여러 수정을 동시에 하지 말고, 원인을 특정한 뒤 하나씩 고칠 것
