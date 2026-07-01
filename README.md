# DiffSBDD in docker-compose

Same shape as the reinvent4-mol2mol setup: a shared `x-diffsbdd-base` anchor, a
`builder` service that produces `htvs-diffsbdd:latest`, runtime knobs in
`host.env`, and a `docker-compose.gpu.yml` override that switches the build to
the CUDA env and reserves the GPU. The difference is that DiffSBDD's entry
points are CLI scripts, so each service reads its flags from `host.env` instead
of a TOML.

## Files

| File | Role |
| --- | --- |
| `Dockerfile` | miniforge base; builds the `diffsbdd` conda env from a selectable env file |
| `docker/environment.gpu.yaml` | upstream CUDA 11.8 env (vendored copy) |
| `docker/environment.cpu.yaml` | CPU-only variant of upstream `environment.yaml` |
| `docker-compose.yml` | base anchor + one service per script (`generate`, `inpaint`, `optimize`, `test`, `train`), plus `fetch-checkpoints` and `shell` |
| `docker-compose.gpu.yml` | GPU override: CUDA env build arg + `nvidia` device reservation |
| `host.env.sample` | thread limits + per-script parameters |
| `.env.gpu.sample` | rename to `.env` to make compose GPU-aware automatically |
| `scripts/fetch_checkpoints.sh` | pull a pretrained checkpoint from Zenodo |

## Install

The **image build** is self-contained — both conda env files are vendored under
`docker/`, so `docker compose build` works from the bundle folder as-is and does
not need DiffSBDD's `environment.yaml` in the context.

To **run** the scripts, the DiffSBDD source (`generate_ligands.py`, `example/`,
etc.) must be present in the workspace.

(If you only want to build and smoke-test the image, you can build from the
bundle folder directly; you just won't have the scripts to run until the source
is mounted.)

## Build

```bash
# CPU (default) — builds on any box, no CUDA needed
docker compose build

# GPU — CUDA 11.8 env + device reservation
docker compose -f docker-compose.yml -f docker-compose.gpu.yml build
# ...or rename .env.gpu.sample to .env and just:  docker compose build
```

The CPU/GPU switch is the `ENV_FILE` build arg (which whole conda env file to
solve), the analogue of `TORCH_INDEX` in the reinvent setup. It changes because
`pytorch-scatter` has to match the exact torch + CUDA build, so swapping a
single pip index isn't enough.

## Get a checkpoint

```bash
docker compose run --rm fetch-checkpoints
```

Downloads `CKPT_URL` to `CKPT` (defaults to the conditional full-atom
CrossDocked model). Point those two vars at any of the eight Zenodo checkpoints
to grab a different one.

## Run

Each service runs one script with parameters from `host.env`:

```bash
docker compose run --rm generate    # generate_ligands.py  (de novo design)
docker compose run --rm inpaint     # inpaint.py           (scaffold/linking)
docker compose run --rm optimize    # optimize.py          (QED/SA evolution)
docker compose run --rm test        # test.py              (whole test set)
docker compose run --rm train       # train.py             (training)
docker compose run --rm shell       # interactive bash in the diffsbdd env
```

`--rm` runs it in the foreground and cleans up — good for the example-sized
jobs. For a long `train` or `test` run, detach and let the daemon own it so an
SSH drop doesn't kill it:

```bash
docker compose up -d train
docker compose logs -f train
docker compose stop train     # graceful SIGTERM, waits up to 180s for a checkpoint
```

Override a single parameter without editing `host.env`:

```bash
N_SAMPLES=100 GENERATE_EXTRA="--sanitize --timesteps 50" \
  docker compose run --rm generate
```

## Notes

- **Device flag.** DiffSBDD's scripts auto-select CUDA when torch sees a GPU and
  fall back to CPU otherwise, so the CPU image needs no extra flag. If you build
  the GPU image but want to force CPU for a quick test, build the CPU image
  instead — there's no per-run device switch in the upstream scripts.
- **`data/` for `test`/`train`.** Those services expect a processed dataset
  (`process_crossdock.py` / `process_bindingmoad.py` output) under the mounted
  workspace. Set `TEST_DIR` / `TRAIN_CONFIG` accordingly.
- **W&B.** `train.py` uses Weights & Biases. Add `WANDB_API_KEY=...` and/or
  `WANDB_MODE=offline` to `host.env` if you train.
- **Image name.** `htvs-diffsbdd:latest`, to sit alongside your
  `htvs-pipeline:latest` and `htvs-redock:latest` images.
