import subprocess
import os

def _get_version():
    """Get version from package metadata or git tags"""
    
    # 1. Try to get from installed package metadata (for installed packages)
    try:
        from importlib.metadata import version
        return "v" + version("ebolaseq")
    except ImportError:
        # Python < 3.8 fallback
        try:
            import pkg_resources
            return "v" + pkg_resources.get_distribution("ebolaseq").version
        except:
            pass
    except:
        pass
    
    # 2. Try to get the latest git tag directly (works in development)
    try:
        result = subprocess.run(['git', 'describe', '--tags', '--abbrev=0'], 
                              capture_output=True, text=True, 
                              cwd=os.path.dirname(__file__))
        if result.returncode == 0:
            return result.stdout.strip()
    except:
        pass
    
    # 3. Fallback: try to get version from git in parent directories
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
    
    # 4. Final fallback if all methods fail
    return "unknown"

__version__ = _get_version()
