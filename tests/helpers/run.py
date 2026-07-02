import subprocess


def run_binary(binary, args, cwd, env, timeout=300, np=None):
    """Run a binary in cwd with env, asserting a clean exit.

    args entries are stringified so Paths and numbers can be passed directly.
    Pass np=<ranks> to launch under `mpirun -np <ranks>`. On failure the
    assertion carries the exit code and the tail of stderr.
    """
    cmd = ["mpirun", "-np", str(np)] if np else []
    cmd += [str(binary), *(str(a) for a in args)]
    proc = subprocess.run(
        cmd,
        cwd=str(cwd),
        env=env,
        capture_output=True,
        text=True,
        timeout=timeout,
    )
    assert proc.returncode == 0, (
        f"{binary} exited {proc.returncode}\n"
        f"args: {args}\n--- stderr (tail) ---\n{proc.stderr[-2000:]}"
    )
    return proc
