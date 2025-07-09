import subprocess
import os

def _get_version():
    """Get version from setuptools_scm generated file or git tags"""
    try:
        # First try the setuptools_scm generated version file (for installed packages)
        from ._version import version as scm_version
        # Extract clean version from setuptools_scm format
        import re
        match = re.match(r'^v?(\d+\.\d+\.\d+)', scm_version)
        if match:
            return f"v{match.group(1)}"
        return scm_version
    except ImportError:
        # Fallback to git tags (for development)
        try:
            result = subprocess.run(['git', 'describe', '--tags', '--abbrev=0'], 
                                  capture_output=True, text=True, 
                                  cwd=os.path.dirname(__file__))
            if result.returncode == 0:
                return result.stdout.strip()
        except:
            pass
        # Final fallback
        return "v0.1.3"

__version__ = _get_version()
