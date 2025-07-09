import subprocess
import os

def _get_version():
    """Get version from git tags - pure git-based approach with fallback"""
    try:
        # Try to get the latest git tag directly (works in development)
        result = subprocess.run(['git', 'describe', '--tags', '--abbrev=0'], 
                              capture_output=True, text=True, 
                              cwd=os.path.dirname(__file__))
        if result.returncode == 0:
            return result.stdout.strip()
    except:
        pass
    
    # Fallback for installed packages: try to get version from git in parent directories
    try:
        # Look for git repo in current or parent directories
        current_dir = os.path.dirname(os.path.abspath(__file__))
        for _ in range(5):  # Check up to 5 parent directories
            if os.path.exists(os.path.join(current_dir, '.git')):
                result = subprocess.run(['git', 'describe', '--tags', '--abbrev=0'], 
                                      capture_output=True, text=True, cwd=current_dir)
                if result.returncode == 0:
                    return result.stdout.strip()
            current_dir = os.path.dirname(current_dir)
    except:
        pass
    
    # Final fallback: if no git available, return a reasonable default
    # This happens when the package is installed and no git repo is available
    return "v0.1.3"  # Current known version

__version__ = _get_version()
