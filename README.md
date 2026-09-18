# code_demo — 재현 코드

문서(00_STATE / 01_PLAN / 02_RESULTS_LOG / HANDOFF)는 이 repo에 없다. project knowledge가 유일한 출처다.
이 repo는 코드만 담는다. 결과 파일(.log / .json / .npy)은 커밋하지 않는다.

```
code_demo/
  code/       스크립트 전부. 이 안에서 실행한다
  results/    실행 로그를 여기로 리다이렉트한다 (git 제외)
```

## 실행

```
pip install -r requirements.txt
cd code
python3 sanity.py                                                           # R6 구현 검증, 약 3분
python3 p26_center_test.py 1.5 300 16                                       # #4 [EXP], 약 12분 CPU
python3 run_ser.py 8 1 4 4000 11              > ../results/ser_8_1_4.log
python3 run_bler.py 8 1 4 60 300 900 5 16,32 1.0:3.0:0.5 > ../results/bler_8_1_4.log
```

R6 스크립트는 `cubesplit_gap.py`를 같은 디렉터리에서 import하므로 `code/` 안에서 실행한다.

## 파일

| 파일 | 항목 |
|---|---|
| code/verify_posterior.py | R1 정식화 검증 |
| code/normalizer_and_sampler.py | R2 정규화 상수, gold-standard MH sampler |
| code/scale_and_cost.py | R3 MACG 스케일, R4 비용 프로파일링 |
| code/gap_experiment.py | R5 스펙트럼 폐기 비용 (대용 shortcut) |
| code/cubesplit_gap.py | R6 라이브러리: Cube-split 원논문 구성·복호기·식 (26), exact ML, LDPC |
| code/sanity.py | R6 구현 검증 |
| code/run_ser.py | R6 SER 스윕. `run_ser.py T B0 N [nblocks] [seed]` |
| code/run_bler.py | R6 BLER 스윕. `run_bler.py T B0 N nblk [ncw_min] [ncw_max] [seed] [etas] [lo:hi:step]` |
| code/p26_center_test.py | R6 (e), 미실행. `p26_center_test.py rho_db ncw [eta]` |

의존성은 numpy, scipy뿐. GPU 불필요.

## 결과 반환

stdout 원문(요약 금지), 스크립트가 남긴 JSON, 실행한 명령 그대로, 코드 수정이 있으면 diff, 환경(CPU/GPU, python·numpy·scipy 버전).
