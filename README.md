![Title](https://github.com/meeravarshneya1234/IKs_stabilizes_APs/blob/master/images/Title.jpg)  

This folder contains source code to reproduce the figures and supplemental material for "IKs stabilizes ventricular action potential" by Varshneya, Devenyi, and Sobie. This work was performed in the Cardiac Systems Pharmacology Lab of Dr. Eric A. Sobie in the Department of Pharmacological Sciences at the Icahn School of Medicine Mount Sinai. 

## Requirements
MATLAB - version 2014 or higher; 2016a was used to run all simulations.

## Usage 
Figures can be recreated by running the "FigureX.m" files. Functions used within each m-file are within the functions folder and have been extensively commented for easy usability. Each model's dydt function was recreated using the equations from the original manuscript. The Heijman model's code was downloaded and edited from the Rudy Lab Website. 

## Ventricular Myocyte Models 
We analyzed ten published and well-known mathematical models describing ventricular myocytes from human, rabbit, canine, and guinea pig. 

| Canine        | Rabbit         | Guinea Pig     | Human          |
| :---:         |     :---:      |          :---: |          :---: |
| [**Fox Model**](https://www.physiology.org/doi/abs/10.1152/ajpheart.00612.2001?url_ver=Z39.88-2003&rfr_id=ori%3Arid%3Acrossref.org&rfr_dat=cr_pub%3Dpubmed)  | [**Shannon Model**](https://www.sciencedirect.com/science/article/pii/S0006349504738023?via%3Dihub)    | [**Livshitz Model**](https://www.sciencedirect.com/science/article/pii/S000634950901159X?via%3Dihub)    |[**TT04 Model**](https://www.physiology.org/doi/abs/10.1152/ajpheart.00794.2003?url_ver=Z39.88-2003&rfr_id=ori%3Arid%3Acrossref.org&rfr_dat=cr_pub%3Dpubmed)      |
| [**Hund Model**](https://www.ahajournals.org/doi/10.1161/01.CIR.0000147231.69595.D3?url_ver=Z39.88-2003&rfr_id=ori%3Arid%3Acrossref.org&rfr_dat=cr_pub%3Dpubmed)     |       | [**Devenyi Model**](https://physoc.onlinelibrary.wiley.com/doi/abs/10.1113/JP273191)      |[**TT06 Model**](https://www.physiology.org/doi/abs/10.1152/ajpheart.00109.2006?url_ver=Z39.88-2003&rfr_id=ori%3Arid%3Acrossref.org&rfr_dat=cr_pub%3Dpubmed)   |
| [**Heijman Model**](https://www.sciencedirect.com/science/article/pii/S002228281100068X?via%3Dihub)     |       |    |[**Grandi Model**](https://www.sciencedirect.com/science/article/pii/S0022282809004295?via%3Dihub)     |
|      |        |  |  [**O'Hara Model**](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002061) |

## Abstract
![Abstract](https://github.com/meeravarshneya1234/IKs_stabilizes_APs/blob/master/images/Graphical_Abstract.jpg)  

The slow and rapid delayed rectifier K+ currents (IKs and IKr, respectively) are responsible for repolarizing the ventricular action potential(AP) and preventing abnormally long APs that may lead to arrhythmias. Although differences in biophysical properties of the two currents have been carefully documented, the respective physiological roles of IKr and IKs are less established. In this study, we sought to understand the individual roles of these currents and quantify how effectively each stabilizes the AP and protects cells against arrhythmias across multiple species. We compared 10 mathematical models describing ventricular myocytes from human, rabbit, canine, and guinea pig. We examined variability within heterogeneous cell populations, tested the susceptibility of cells to proarrhythmic behavior, and studied how IKs and IKr responded to changes in the AP. We found that: (1) models with higher baseline IKs exhibited less cell-to-cell variability in action potential duration (APD); (2) models with higher baseline IKs were less susceptible to early afterdepolarizations (EADs) induced by depolarizing perturbations; (3) as APD is lengthened, IKs increases more profoundly than IKr, thereby providing negative feedback that resists excessive AP prolongation; and (4) the increase in IKs that occurs during β-adrenergic stimulation is critical for protecting cardiac myocytes from EADs under these conditions. Slow delayed rectifier current is uniformly protective across a variety of cell types. These results suggest that IKs enhancement could potentially be an effective antiarrhythmic strategy.

<p align="center">
  <b>"IKs provides better negative feedback than IKr"</b><br>
</p>

![GIF](https://github.com/meeravarshneya1234/IKs_stabilizes_APs/blob/master/images/IKsFeedbackLoop.gif)  

## Questions/Contributions
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

## Citation 
M. Varshneya, R. A. Devenyi, and E. A. Sobie, “Slow Delayed Rectifier Current Protects Ventricular Myocytes From Arrhythmic Dynamics Across Multiple Species,” Circ. Arrhythmia Electrophysiol., vol. 11, no. 10, Oct. 2018.
