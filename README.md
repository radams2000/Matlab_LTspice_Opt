# Matlab_LTspice_Opt
Matlab optimizer with LTspice

Written by Bob Adams, Analog Devices Inc.
## LTSpice_opt is a Matlab program that uses an iterative optimization approach to design analog filters. 
It is designed to be used in conjunction with the popular circuit simulator LTspice. It embeds an LTSpice simulation inside the powerful Matlab nonlinear least-squares optimizer (lsqnonlin). 
It works as follows;

1) The user provides a target frequency and/or phase response in Matlab, and a circuit topology in LTspice with some initial component values.
  
2) The user provides a list of which circuit instances the optimizer is allowed to vary.
  
3) The optimizer then iteratively adjusts those component values, running a simulation for every pass through the lsqnonlin algorithm, in an attempt to reduce the error between the target frequency (and/or phase) response and the simulated response.
  
4) Once the optimizer has finished, a new schematic is generated with the optimized component values. During the schematic generation process, each component value is quantized to a user-defined tolerance.

Why is this capability useful? Don't we already know how to design filters?

Traditional filter design uses standard circuit topologies such as Sallen-and-Key or multiple-feedback op-amp based active filters. In these cases, given a "standard" filter shape such a Butterworth or Chebychev, a high-order filter may be factored into 2nd-order sections, and an op-amp circuit can be used for each of those sections. This design procedure is quite straightforward and has not changed for many years.


However, there are many cases where this approach fails;

***The desired filter shape is not a traditional shape such as Butterworth or Chebychev***.

 For example, the filter may need to compensate for some other part of the system that has a non-flat frequency response, while simultaneously attenuating other frequency regions.

 ***The need to reduce power/area by combining multiple filter sections into a single op-amp circuit***.
 
This leads to non-conventional circuit topologies that have very messy closed-loop formulas, and it becomes very difficult to solve for the component values. 

 ***The application operates at frequencies where finite op-amp gain-bandwidth degrades the frequency response***.
 
Calculating the effects of finite gain-bandwidth on the frequency response is quite complicated, especially in cases where the gain/phase response deviates from the traditional single-pole model. By running optimizer simulations in LTSpice, the actual target op-amp may be included in the simulation. This yields a solution that inherently attempts to compensate for finite gain-bandwidth effects.

The software may be downloaded using the link here.

Download the zip file, and unzip in a directory on your computer. Move all the Matlab ".m" files to a directory where you will run Matlab. Move the LTspice ".asc" files to a directory where you will store your schenatics. See the Readme.pdf file for detailed instruction on how to run the Matlab program. 


