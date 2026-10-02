import curses
import os
from logfile_manager import logfile_manager

def init_curses():
    stdscr = curses.initscr()
    curses.savetty()
    curses.start_color()
    curses.use_default_colors()
    curses.init_pair(1, curses.COLOR_WHITE,  curses.COLOR_BLACK)
    curses.init_pair(2, curses.COLOR_YELLOW, curses.COLOR_BLACK)
    curses.init_pair(3, curses.COLOR_CYAN,   curses.COLOR_BLACK)
    curses.init_pair(4, curses.COLOR_RED,    curses.COLOR_BLACK)
    curses.init_pair(5, curses.COLOR_GREEN,  curses.COLOR_BLACK)
    curses.init_pair(6, curses.COLOR_BLUE,   curses.COLOR_BLACK)
    curses.curs_set(0)
    curses.noecho()
    curses.cbreak()
    stdscr.keypad(True)
    return stdscr

def close_curses(stdscr):
    curses.nocbreak()
    stdscr.keypad(False)
    curses.echo()    
    curses.curs_set(1)
    curses.resetty()
    curses.endwin()

COLOR_MAP = {
    'd': 0,  # Default
    'w': 1,  # White
    'y': 2,  # Yellow
    'c': 3,  # Cyan
    'r': 4,  # Red
    'g': 5,  # Green
    'b': 6   # Blue
}

def log(panel, indent, text, color='c'):
    color_pair = COLOR_MAP.get(color, 0)
    height, width = panel.getmaxyx()
    wrapped_lines = wrap_text_with_indent(text, width, indent)
    logfile = logfile_manager.logfile
    with open(logfile, 'a') if logfile else open(os.devnull, 'w') as log_file:
        for line in wrapped_lines:
            try:
                panel.addstr(panel.getyx()[0], 0, f"{line}\n", curses.color_pair(color_pair))
                if logfile:
                    log_file.write(f"{line}\n")
            except curses.error:
                pass
    panel.refresh()

def wrap_text_with_indent(text, width, indent):
    words = text.split(' ')
    wrapped_lines = []
    current_line = ' ' * indent

    for word in words:
        if len(current_line) + len(word) + 1 > width:
            wrapped_lines.append(current_line)
            current_line = ' ' * indent + word + ' '
        else:
            current_line += word + ' '

    wrapped_lines.append(current_line)
    return wrapped_lines
