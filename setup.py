from setuptools import setup


def branch_aware_local(version):
    """Include branch name in local version segment for non-master branches.

    Examples:
      Tagged commit (any branch):      1.2.0
      master, untagged:                1.3.0.dev79+gabcdef1
      dev branch, untagged:            1.3.0.dev79+dev.gabcdef1
      feature/foo branch, untagged:    1.3.0.dev79+feature.foo.gabcdef1
    """
    if version.exact:
        return ""
    node = f"g{version.node[:6]}" if version.node else ""
    branch = version.branch
    if branch and branch not in ("master", "main", "HEAD"):
        safe = branch.replace("/", ".").replace("-", ".").replace("_", ".")
        return f"{safe}.{node}" if node else safe
    return node


setup(
    use_scm_version={
        "write_to": "easydock/_version.py",
        "version_scheme": "guess-next-dev",
        "local_scheme": branch_aware_local,
        "fallback_version": "0.0.0+unknown",
    }
)
