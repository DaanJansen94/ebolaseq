import subprocess

def _get_version_from_git():
    """Get the latest git tag version"""
    result = subprocess.run(['git', 'describe', '--tags', '--abbrev=0'], 
                          capture_output=True, text=True, cwd=__file__.rsplit('/', 1)[0])
    if result.returncode == 0:
        return result.stdout.strip()
    raise RuntimeError("Could not determine version from git tags")

__version__ = _get_version_from_git() 