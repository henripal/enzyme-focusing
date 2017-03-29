# Trajectory analysis

The experimental data is a set of videos with different experimental
conditions. C=0 is the green channel, imaging HK. C=1 is the red channel,
imaging Aldolase.

For each video, we identify the Hexokinase trajectories. Then, for each
hexokinase trajectory, we plot all nearby Aldolase trajectory and plot the
result.

To do: calculate a significance test that shows the trajectories are related
(although it is obvious from visual inspection).

## Recap of the algorithmic procedure

The steps are shown in [this notebook](../rev_notebooks/particle_expl.ipynb).

- locate step: the image goes through a band pass filter, a threshold filter,
  then the peaks are located and fit to a gaussian (using trackpy)
- trajectory step: the features are then linked from frame to frame into
  trajectories.

The parameters of the above are manually tuned to produce meanigful looking
trajectories, then are fixed for all videos, except the 20X video. (to do:
store the parameters!)

The *blue* curves are the hexokinase trajectories, the *red* curves are the
aldolase trajectories.

### Glucose 1 video results

![](../img/gl_1.png)

### Glucose 2 video results

![](../img/gl_2.png)

### Glucose 3 video results

![](../img/gl_3.png)

### Glucose 4 video results

![](../img/gl_4.png)

### Glucose 5 video results

![](../img/gl_5.png)

### Control results

The same procedure applied to the control videos yield either one or zero
hexokinase trajectories, with no overlapping or nearby aldolase trajectories,
ending the procedure.

For example here are the results for the l-glucose videos:
![](../img/lglu_001.png)
![](../img/lgl.png)
