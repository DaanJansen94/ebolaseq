import subprocess
import os

def _get_version():
    """Get version from git tags - always show current tag, not next version"""
    try:
        # First try to get the latest git tag directly (for both dev and installed)
        result = subprocess.run(['git', 'describe', '--tags', '--abbrev=0'], 
                              capture_output=True, text=True, 
                              cwd=os.path.dirname(__file__))
        if result.returncode == 0:
            return result.stdout.strip()
    except:
        pass
    
    try:
        # Fallback: try setuptools_scm but extract the actual current version
        from ._version import version as scm_version
        import re
        
        # If it's a dev version (e.g., "0.1.4.dev18+..."), we're still on the previous tag
        dev_match = re.match(r'^v?(\d+)\.(\d+)\.(\d+)\.dev', scm_version)
        if dev_match:
            major, minor, patch = dev_match.groups()
            # Decrement patch version to get the actual current tag
            current_patch = max(0, int(patch) - 1)
            return f"v{major}.{minor}.{current_patch}"
        
        # If it's a clean version, use it as-is
        clean_match = re.match(r'^v?(\d+\.\d+\.\d+)$', scm_version)
        if clean_match:
            return f"v{clean_match.group(1)}"
            
        return scm_version
    except ImportError:
        pass
    
    # If we can't determine version, something is wrong
    raise RuntimeError("Could not determine version from git tags or setuptools_scm")

__version__ = _get_version()
