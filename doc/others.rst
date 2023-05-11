.. _other:

Additionnal Functionalities
===========================

Here are more functionalities implemented in the tool.

Inputs analysis Tool
--------------------

For computing a leaf area distribution, you can run s5 or s2v tools, which was developped in previous projects. 
This tools are accessible with the methods :meth:`s5` and :meth:`s2v`.
Those tools will create output files with leaf angle distribution, leaf area density and 
other informations related to turbid medium, from a triangles mesh inside a LightVegeManager object.


Visualisation
-------------

VTK files
*********

You can create VTK files from the geometric information you have in the inputs. Files can be either triangles or voxels mesh. Available methods
    
    * :meth:`VTK_nolight`: write only geometric informations from the common scene
    
    * :meth:`VTK_light`: write geometric informations associated with radiation values for each element
    
    * :meth:`VTK_sun`: write a line representing the sun position

PlantGL visualizer
******************

Currently in development.