If you uncomment the line 
if (iter%iter_number == 0) bov_window_screenshot(window, screenshot_name);
in the routine 
print_particles > display_particles
the code will output you "*.ppm" images every "iter_number" iterations.

To make a movie from those images, use the 3 scripts in this folder in the order indicated by the numbers:
1. Convert the strange "*.ppm" format into "*.jpg"
2. Rename correctly the "*.jpg" format such that the images are assembled in the right order to make the movie at the next step
3. Convert the sequence of "*.jpg" images into a "*.mpg" movie (rem: need to convert that to "*.mp4" format to be uploaded in a PowerPoint presentation)
