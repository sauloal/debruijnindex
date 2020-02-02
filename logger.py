import sys

LOG_LEVEL_DEBUG   = 0
LOG_LEVEL_LOG     = 2
LOG_LEVEL_INFO    = 4
LOG_LEVEL_WARNING = 6
LOG_LEVEL_ERROR   = 8

LOG_LEVEL_DEFAULT = LOG_LEVEL_INFO
LOG_LEVEL_DEFAULT_VERBOSITY = LOG_LEVEL_INFO

LOG_LEVEL = [
    "DEBUG  ", # 0
    "DEBUG  ", # 1
    "LOG    ", # 2
    "LOG    ", # 3
    "INFO   ", # 4
    "INFO   ", # 5
    "WARNING", # 6
    "WARNING", # 7
    "ERROR  ", # 8
    "ERROR  ", # 9
]

LOG_COLORS = [
    '\033[95m', # 0 DEBUG
    '\033[95m', # 1 DEBUG

    '\033[94m', # 2 LOG
    '\033[94m', # 3 LOG
    
    '\033[92m', # 4 INFO
    '\033[92m', # 5 INFO

    '\033[93m', # 6 WARNING
    '\033[93m', # 7 WARNING
    
    '\033[91m', # 8 ERROR
    '\033[91m', # 9 ERROR
]

LOG_COLOR_ENDC      = '\033[0m'
LOG_COLOR_BOLD      = '\033[1m'
LOG_COLOR_UNDERLINE = '\033[4m'

sysprint = print

def print(*args, **kwargs):
    level = LOG_LEVEL_DEFAULT
    verbosity = LOG_LEVEL_DEFAULT_VERBOSITY
    
    if "level" in kwargs:
        level = kwargs["level"]
        del kwargs["level"]

    if "verbosity" in kwargs:
        verbosity = kwargs["verbosity"]
        del kwargs["verbosity"]

    if level >= 6:
        if "file" not in kwargs:
            kwargs["file"] = sys.stderr

    if verbosity <= level:
        msg = [LOG_COLORS[level], LOG_LEVEL[level], ": ", LOG_COLOR_ENDC] + list(args) + [LOG_COLOR_ENDC]
        sysprint(*msg, **kwargs)


def print_debug(*args, **kwargs):
    kwargs["level"] = LOG_LEVEL_DEBUG
    print(*args, **kwargs)
    
def print_log(*args, **kwargs):
    kwargs["level"] = LOG_LEVEL_LOG
    print(*args, **kwargs)

def print_info(*args, **kwargs):
    kwargs["level"] = LOG_LEVEL_INFO
    print(*args, **kwargs)

def print_warning(*args, **kwargs):
    kwargs["level"] = LOG_LEVEL_WARNING
    print(*args, **kwargs)

def print_error(*args, **kwargs):
    kwargs["level"] = LOG_LEVEL_ERROR
    print(*args, **kwargs)
    