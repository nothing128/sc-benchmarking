import sys, os, gc
import polars as pl
from contextlib import contextmanager
from timeit import default_timer
sys.path.append('projects/def-wainberg/karbabi/utils')
from single_cell import SingleCell

# region functions

class TimerCollection:
    def __init__(self, silent=True):
        self.timings = {}
        self.silent = silent
        
    def __call__(self, message):
        start = default_timer()
        @contextmanager
        def timer():
            if not self.silent:
                print(f'{message}...')
            try:
                yield
                aborted = False
            except Exception as e:
                aborted = True
                raise e
            finally:
                self.timings[message] = {
                    'duration': default_timer() - start,
                    'aborted': aborted
                }
        return timer()
    
    def print_summary(self, sort=True):
        print('\n--- Timing Summary ---')
        if sort:
            timings_items = sorted(
                self.timings.items(), 
                key=lambda x: x[1]['duration'], 
                reverse=True)
        else:
            timings_items = list(self.timings.items())
        total_time = sum(info['duration'] for _, info in timings_items)
        
        for message, info in timings_items:
            duration = info['duration']
            percentage = (duration / total_time) * 100 if total_time > 0 else 0
            status = 'aborted after' if info['aborted'] else 'took'
            time_str = self._format_time(duration)
            print(f'{message} {status} {time_str} ({percentage:.1f}%)')
        print(f'\nTotal time: {self._format_time(total_time)}')
    
    def _format_time(self, duration):
        units = [
            (86400, 'd'),
            (3600, 'h'),
            (60, 'm'),
            (1, 's'),
            (0.001, 'ms'),
            (0.000001, 'µs'),
            (0.000000001, 'ns')
        ]
        parts = []
        for threshold, suffix in units:
            if duration >= threshold or \
                (not parts and threshold == 0.000000001):
                if threshold >= 1:
                    value = int((duration // threshold))
                    duration %= threshold
                else:
                    value = int((duration / threshold) % 1000)
                if value > 0 or (not parts and threshold == 0.000000001):
                    parts.append(f'{value}{suffix}')
                if len(parts) == 2:
                    break
        return ' '.join(parts) if parts else 'less than 1ns'
    
# endregion

# endregion

timers = TimerCollection()

with timers('task 1'):

