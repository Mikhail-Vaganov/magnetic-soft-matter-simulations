#Magnetic soft matter simulations

The project contains simulations of a soft matter with magnetic particles.
Magnetically *hard* particles simulated by the Stoner-Wohlfarth model.
Magnetically *soft* particles obey the Froehlich-Kennelly law.

The code is written in Matlab.

### Stoner-Wolfarth model simulation

The Stoner-Wolfarth model of the magnetic hysteresis is implemented in SwParticle class. 
The constructor requires one argument: the angle between the applied field and the easy axis (uniaxial anisotropy) 

```matalb
p = SwParticle(pi/3);
p.Draw('Folder for image representations');
```

### First order reversal curve (FORC) simulation

FORC folder contains simulation of the first order reversal curve experiment.

Simulation of the FORC experiment consists of four steps:

1. FORC processor initialization
2. Preparation of the matter (processing of a hysteresis loop for acceleration of the simulation)
3. FORCs measurement
4. First order reversal curve diagrams calculation
5. Visualizing of the FORC diagrams

The FORC diagrams analysis is implemented by the FORC processors, for instance, by objects of PikeFORC class.
FORC analysis can be applied to a matter simulating object which implements iMatter interface. 
For example, in order to simulate a single particle model use SingleParticleMatter class, whose constructor takes in turn an implementation of iMagneticParticle interface.

```matalb
p = SwParticle(pi/3);
matter = SingleParticleMatter(p);
forc = PikeFORC(4, -4, 4, matter, resultsFolder);
forc = forc.PrepareMatter();
forc = forc.MagnetizationFORC();
forc = forc.CalculateFORCDistribution();
forc.DrawResults();
```

### Preisach hysteresis model simulation

In the Preisach model of hysteresis magnetic particles represented by the abstracted elements called hysterons, which possess rectangular hysteresis loop.
This conception is implemented by Hysteresis class. 
The constructor requires the switching field values, at which magnetization of the particle jumps to the upper and bottom branches correspondingly.

```matlab
h = Hysteron(4, 0);
h.Draw('Folder for image representations');
```

The model of hysteresis as sum of hysterons is presented in the corresponding Preisach folder

### Usage of this code

I will be glad to collaborate with you on other magnetic sumulations.

If you decide to use parts of my code in your research, please, let me know via account on GitHub or via e-mail sent to mikhail.vaganov.sci@gmail.com

Or, just cite this [article](http://www.sciencedirect.com/science/article/pii/S0304885316319552):
* Vaganov, Linke, Odenbach, Raikher, JMMM 2015


##References
1. Pike, Roberts, Verosub, JAP 1999
2. Winklhofer, Zimanyi, JAP 2006
3. Harrison, Feinberg, GGGE 2008 
