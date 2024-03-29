.. This is A COPY OF the main index.rst file which is rendered into the landing page of your documentation.
   Follow the inline instructions to configure this for YOUR next project.



Welcome to Telescope Scheduler's documentation !
=========================================================
|

A scheduler for observations on a high-cadence telescope.

Main objective:
    - given a list of positions in the sky
    - maximize the number of positions which are observed at least three times
    - respecting a minimum interval between consecutive observations of the same position

The source code is available `here <https://github.com/joaopedroso/telescope_scheduler>`_.

Workflow for adapting the code to a different setting:

1. Specify telescope moving speed in `motion_time.py` and and location, and exposure times in `constants.py`.
2. Implement functions `mk_obs_time()` and `mk_obs_set()`, returning respectively the observation times and the set of sky positions to observe, in file `instance.py`.
3. For real-time updates of positions that cannot be observed (e.g., due to being hidden by clouds), implement function `hidden(sky)` in file `hidden.py`.
4. Adapt file `simulator.py` for your situation.
   

|

.. maxdepth = 1 means the Table of Contents will only links to the separate pages of the documentation.
   Increasing this number will result in deeper links to subtitles etc.

.. Below is the main Table Of Contents
   You have below a "dummy" file, that holds a template for a class.
   To add pages to your documentation:
        * Make a file_name.rst file that follows one of the templates in this project
        * Add its name here to this TOC


.. toctree::
   :maxdepth: 1


.. autosummary:: 
   :toctree: generated

   simulator
   scheduler
   hidden
   instance
   motion_time
   obs_data
   solution
   astro_tools
   constants


.. Delete this line until the * to generate index for your project: * :ref:`genindex`


|

This documentation was last updated on |today|.

.. Finished personalizing all the relevant details? Great! Now make this your main index.rst,
   And run `make clean html` from your documentation folder :)
