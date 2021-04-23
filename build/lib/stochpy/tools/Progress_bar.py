 #! /usr/bin/env python
"""
Progress bar
===========

Written by Maxim Moinat and TR Maarleveld, Amsterdam, The Netherlands
E-mail: tmd200@users.sourceforge.net
Last Change: February 18, 2015
"""

from __future__ import division, print_function, absolute_import
import sys,time

class Progress_bar():
    """
    A progress bar for slow iterative processes.
    Give problems if print statements are in the iterative process, the bar is then not refreshed but recreated
    """
    def __init__(self, cycles_total, total_signs = 20, done_msg = "Info: Done"):      
        self.i = 0        
        self.end = cycles_total
        self.total_signs = total_signs
        
        self.basic_bar = "[%%-%ss] %%d%%%%, Runtime: %%.2f sec" % self.total_signs
        self.done_msg = done_msg
        
        self.t1 = time.time()
        #self.update(add = 0) # Shows bar before start of iteration (no step is added)
        #Without above line, there can be a slight delay with start script and showing of bar. -> Line disabled for compatibility with StochSim
    
    def update(self, add = 1,quiet=False):        
        self.i += add
        #Remove bar if end of loop reached. Note: 100% is never printed, the progressbar immediately ends when reached.
        if self.i >= self.end: # With one cycle the progressbar will not be printed at all, (unless initialized) -> Not printed for compatibility with StochSim
            self.done_bar(quiet)
            return None # Escapes the function
            
        self.frac = self.i/float(self.end)
        self.n_signs = int(self.total_signs * self.frac)        
        self.print_()
        
    def print_(self):        
        self.t2 = time.time()
        self.new_bar = self.basic_bar %('|'*self.n_signs, self.frac*100, self.t2-self.t1)
        sys.stdout.write('\r') # Resets cursor to begin of sentence, which is to be rewritten
        sys.stdout.write(self.new_bar)
        sys.stdout.flush()
        
    def done_bar(self,quiet):
        if self.done_msg == 'time': # Special end message: the simulation time
            self.end_time = time.time()
            self.done_msg = "Info: Simulation time: {0}".format(self.end_time - self.t1)
        
        sys.stdout.write('\r')        
        if not quiet:        
            sys.stdout.write(self.done_msg + ' '*(len(self.basic_bar)+self.total_signs)) # + '\n') # Done message, some excess whitespace to make sure to delete the completed progress bar, and newline to print possible new text on the next line       
        else:     
            sys.stdout.write(' '*(len(self.basic_bar)+self.total_signs)+'\r')         
        sys.stdout.flush()
        
def main():
    NUMBER_OF_EVENTS = 10
    t1 = time.time()
    
    bar = Progress_bar(cycles_total = NUMBER_OF_EVENTS, total_signs = 20, done_msg = 'time')
    for i in range(NUMBER_OF_EVENTS):
        for g in range(1000):
            a = [g]*g**2
            s = 0
            for j in a:
                s+= j
                bar.update()
        
    t2 = time.time()
    
    print(t2 - t1)
    
if __name__ == '__main__':
    main()    
