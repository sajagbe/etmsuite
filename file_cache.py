#!/usr/bin/env python3
"""
File Caching System for ETM Suite
=================================

High-performance caching system to eliminate redundant file I/O operations.
This module provides thread-safe caching for parsed quantum chemistry output files.

Instead of reading the same file 10+ times
for different properties, we read it once and cache the parsed data.
"""

import os
import threading
from functools import lru_cache
from typing import Dict, Any, Optional

try:
    from cclib.io import ccread
    CCLIB_AVAILABLE = True
except ImportError:
    CCLIB_AVAILABLE = False

class FileCache:
    """Thread-safe cache for parsed quantum chemistry files"""
    
    def __init__(self, max_size: int = 128):
        self._cache: Dict[str, Any] = {}
        self._lock = threading.RLock()
        self._max_size = max_size
        self._access_order = []
    
    def get_parsed_file(self, filepath: str) -> Optional[Any]:
        """Get cached parsed file data or parse and cache it"""
        if not os.path.exists(filepath):
            return None
            
        # Use file modification time + path as cache key for integrity
        mtime = os.path.getmtime(filepath)
        cache_key = f"{filepath}:{mtime}"
        
        with self._lock:
            if cache_key in self._cache:
                # Move to end for LRU tracking
                self._access_order.remove(cache_key)
                self._access_order.append(cache_key)
                return self._cache[cache_key]
            
            # Parse file and cache result
            try:
                if CCLIB_AVAILABLE:
                    parsed_data = ccread(filepath)
                else:
                    # For non-cclib fallback, cache the file content
                    with open(filepath, 'r') as f:
                        parsed_data = f.read()
                
                # Implement LRU eviction
                if len(self._cache) >= self._max_size:
                    oldest_key = self._access_order.pop(0)
                    del self._cache[oldest_key]
                
                self._cache[cache_key] = parsed_data
                self._access_order.append(cache_key)
                return parsed_data
                
            except Exception as e:
                print(f"Error parsing file {filepath}: {e}")
                return None
    
    def clear(self):
        """Clear all cached data"""
        with self._lock:
            self._cache.clear()
            self._access_order.clear()
    
    def get_cache_stats(self) -> Dict[str, int]:
        """Get cache statistics"""
        with self._lock:
            return {
                'size': len(self._cache),
                'max_size': self._max_size,
                'hit_ratio': getattr(self, '_hits', 0) / max(getattr(self, '_requests', 1), 1)
            }

# Global cache instance
_global_cache = FileCache()

def get_cached_file(filepath: str) -> Optional[Any]:
    """Get cached parsed file - global convenience function"""
    return _global_cache.get_parsed_file(filepath)

def clear_cache():
    """Clear global cache"""
    _global_cache.clear()

def get_cache_stats() -> Dict[str, int]:
    """Get global cache statistics"""
    return _global_cache.get_cache_stats()

@lru_cache(maxsize=64)
def get_file_content(filepath: str, mtime: float) -> str:
    """LRU cached file content reader (for regex parsing)"""
    try:
        with open(filepath, 'r') as f:
            return f.read()
    except Exception as e:
        print(f"Error reading file {filepath}: {e}")
        return ""

def get_cached_content(filepath: str) -> str:
    """Get cached file content with modification time check"""
    if not os.path.exists(filepath):
        return ""
    
    mtime = os.path.getmtime(filepath)
    return get_file_content(filepath, mtime)
