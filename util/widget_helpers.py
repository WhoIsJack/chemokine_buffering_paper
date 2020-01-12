# -*- coding: utf-8 -*-
"""
Created on Wed Sep 19 19:27:45 2018

@author:    Jonas Hartmann @ Gilmour group @ EMBL Heidelberg

@descript:  Convenient little helpers for ipywidgets hacking.
"""

#------------------------------------------------------------------------------

# IMPORTS

from __future__ import division
import os
import ipywidgets as widgets
import matplotlib.pyplot as plt
from IPython.display import display, clear_output


#------------------------------------------------------------------------------

# Decorator to add a "save figure" button to matplotlib plots,
# including plots created in interactive widgets using ipywidgets.interact!

def savebutton(func):
    """Decorator that adds a "save figure" button to a function that generates
    a matplotlib plot, including functions made into interactive widgets using
    ipywidgets.interact.
    
    NOTE: make sure you don't have `plt.show()` in the figure generation
          function, as that will clear the current figure!
          
    WARNING: This does not properly forward defaults for sliders and such set
             in the function definition. Instead, the defaults are those that
             are automatically determined by ipywidgets.interact!
          
    Examples
    --------
    
    # Without `interact` widget
    @savebutton
    def make_plot():
        t = np.linspace(0, 10, 500)
        y = np.sin(2.0*t)
        plt.plot(t, y, color='red')

    # With `interact` widget
    @interact(freq=(0.5, 10, 0.5),
              color=['red', 'blue', 'green'])
    @savebutton
    def make_plot(freq=2.0, color='red'):
        t = np.linspace(0, 10, 500)
        y = np.sin(freq*t)
        plt.plot(t, y, color=color)
    
    """
        
    def wrapper(**kwargs):
        
        # Prepare textbox (for filename) and button
        textbox = widgets.Text(value='', placeholder='Enter Filename',
                               description='Filename:', disabled=False)
        button  = widgets.Button(description='Save figure!')
        box     = widgets.HBox([textbox, button])
        
        # Callback to save figure when button is clicked
        # TODO: This current uses/creates a specific folder ('figures') for
        #       saving. It would be great to make this more general, perhaps 
        #       by triggering a standard file saving dialogue to open!
        def on_button_clicked(b):
            if textbox.value:
                if not os.path.isdir('figures'):
                    os.mkdir('figures')
                figpath = str(textbox.value)
                if not figpath.endswith('.pdf'):
                    figpath += '.pdf'
                figpath = os.path.join('figures', figpath)
                b.fig.savefig(figpath, bbox_inches='tight', transparent=True)
                print "Saved figure as '%s'" % figpath
                textbox.value=''
        button.on_click(on_button_clicked)
        
        # Run wrapped function to generate figure
        func(**kwargs)
        
        # Update figure in button
        button.fig = plt.gcf()
        
        # Display textbox and button
        display(box)
        
    # Done!
    return wrapper  


#------------------------------------------------------------------------------



