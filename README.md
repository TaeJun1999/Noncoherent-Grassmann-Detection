# Noncoherent-Grassmann-Detection

Rayleigh block-fading noncoherent MIMO의 Grassmann posterior를 생성 추론으로 푸는 연구.
코드, 원시 로그, 연구 문서의 **원본**이다. Claude project와 GPT project의 문서는 이곳의 사본이다.

## 구조

```
code/                  스크립트 전부. R6 스크립트는 cubesplit_gap.py를 같은 디렉터리에서
                       import하므로 code/ 안에서 실행한다
docs/                  연구 문서
results/<ID>_<설정>/   실행 원시 로그 (cmd.txt, env.txt, sanity.log, stdout.log)
```

## 문서

| 파일 | 역할 |
|---|---|
| `docs/00_STATE.md` | 현재 상태. 다른 문서와 충돌하면 이 파일이 이긴다 |
| `docs/01_PLAN.md` | 정식화·주장·선행연구·baseline·설계. `⚠️ STALE` 절은 인용하지 않는다 |
| `docs/02_RESULTS_LOG.md` | 실행한 실험의 기록. **인용 가능한 숫자의 유일한 출처** |
| `docs/03_LIT_LOG.md` | 검색으로 확인한 선행연구. **논문 인용의 유일한 출처** |
| `docs/derivations/D<n>_*.md` | GPT project 산출물. 제안이며 검산 전에는 확정이 아니다 |
| `docs/handoff/*.md` | chat 간 인수인계 |
| `docs/sessions/*.txt` | 서버 실행 세션 트랜스크립트 |
| `docs/EXPERIMENTS.md` | 실행 대장. 한 실행 = 한 행. 인용 출처가 아니다 |

문서 4개는 마일스톤마다 **통째로** 교체한다. 부분 수정하지 않는다.

## 실행

파이썬은 절대경로로 부른다 (`CLAUDE.md`). PATH의 `python3`에는 numpy가 없다.

```
PY=~/miniforge3/envs/torch/bin/python
cd code

$PY sanity.py                                                            # R6 구현 검증
$PY run_ser.py 8 1 4 4000 11                       > ../results/ser_8_1_4.log
$PY run_bler.py 8 1 4 60 300 900 5 16,32 1.0:3.0:0.5 > ../results/bler_8_1_4.log
$PY p26_center_test.py 1.5 300 16                                        # R6 (e)
```

의존성은 numpy와 scipy뿐이다. R6 재현은 1코어에서 약 50분.

## 기록 규칙

- 숫자는 로그 **원문**에서만 옮긴다. 요약본에서 옮기지 않는다
- 결과를 커밋할 때 `code_commit`과 `result_commit`을 함께 남긴다
- `lists_*.npy`, 체크포인트, `logs/`는 커밋하지 않는다 (`.gitignore`)
