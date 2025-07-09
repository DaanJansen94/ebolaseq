import subprocess
import os

def _get_version():
    """Get version directly from git tags only"""
    
    # Try to get the latest git tag directly (works in development)
    try:
        result = subprocess.run(['git', 'describe', '--tags', '--abbrev=0'], 
                              capture_output=True, text=True, 
                              cwd=os.path.dirname(__file__))
        if result.returncode == 0:
            return result.stdout.strip()
    except:
        pass
    
    # Fallback: try to get version from git in parent directories
    try:
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
    
    # No fallback - if git is not available, return unknown
    return "unknown"

__version__ = _get_version()
