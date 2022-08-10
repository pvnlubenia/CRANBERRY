========================================================

  CRANBERRY: Chemical ReAction Network numBERs summaRY

========================================================

Matlab was used to develop the function used here.
The function networkNumbers returns the values of the following network numbers of a chemical reaction network:

     - Species (m)
     - Complexes (n)
     - Reactant complexes (n_r)
     - Reversible reactions (r_rev)
     - Irreversible reactions (r_irrev)
     - Reactions (r) = r_irrev + 2 r_rev
     - Linkage classes (l)
     - Strong linkage classes (sl)
     - Terminal linkage classes (t)
     - Rank (s)
     - Reactant rank (q)
     - Deficiency (delta) = n - l - s
     - Reactant deficiency (delta_p) = n_r - q

Furthermore, the output variable 'model' allows the user to view the complete network with all the species listed in the 'species' field of the structure 'model'.

Parts of the code were patterned from the file model_analysis.m which is part of the ERNEST toolbox for chemical chemical reaction network theory [3]. The toolbox can be downloaded from https://www.sissa.it/fa/altafini/papers/SoAl09/.

The codes for the computation for the number of reversible reaction, irreversible reactions, reactions, reactant rank, and reactant deficiency are unique from the author. Discussions on reactant rank and reactant deficiency can be found on [2].



=================================
How to fill out 'model' structure
=================================

'model' is the input for the function networkNumbers. It is a structure, representing the CRN, with the following fields:

   - id: name of the model
   - species: a list of all species in the network; this is left blank since incorporated into the function is a step which compiles all species used in the model
   - reaction: a list of all reactions in the network, each with the following subfields:
        - id: a string representing the reaction
        - reactant: has the following further subfields:
             - species: a list of strings representing the species in the reactant complex
             - stoichiometry: a list of integers representing the stoichiometry of each species in the reactant complex (listed in the same order of the species)
        - product: has the following further subfields:
             - species: a list of strings representing the species in the product complex
             - stoichiometry: a list of integers representing the stoichiometry of each species in the product complex (listed in the same order of the species)
        - reversible: has the value 'true' or 'false' indicating if the reaction is reversible or not, respectively

To fill out the 'model' structure, write a string for 'model.id': this is just to put a name to the network. To add the reactions to the network, use the function addReaction where the output is 'model'. addReaction is developed to make the input of reactions of the CRN easier than the input in [3]:

   addReaction
      - OUTPUT: Returns a structure called 'model' with added field 'reaction' with subfields 'id', 'reactant', 'product', 'reversible', and 'kinetic'. The output variable 'model' allows the user to view the network with the added reaction.
      - INPUTS
           - model: a structure, representing the CRN
           - id: visual representation of the reaction, e.g., reactant -> product (string)
           - reactant_species: species of the reactant complex (cell)
           - reactant_stoichiometry: stoichiometry of the species of the reactant complex (cell)
           - reactant_kinetic: kinetic orders of the species of the reactant complex (array)
           - product_species: species of the product complex (cell)
           - product_stoichiometry: stoichiometry of the species of the product complex (cell)
           - product_kinetic: "kinetic orders" of the species of the product complex, if the reaction is reversible (array); if the reaction in NOT reversible, leave blank
           - reversible: logical; whether the reaction is reversible or not (true or false)
      * Make sure the function addReaction is in the same folder/path being used as the current working directory.


========
Examples
========

2 examples are included in this folder, all based on [1]:

   - Example 1: ARL3-S Molecular network of Leishmaniasis infectious disease

   - Example 2: ECJ3-G Yeast trehalose system



===================
Contact Information
===================

For questions, comments, and suggestions, feel free to contact me at pvnlubenia@yahoo.co.uk.


- Patrick Lubenia (18 July 2022)



==========
References
==========

   [1] Arceo C, Jose E, Lao A, Mendoza E (2016) Reaction networks and kinetics of biochemical systems (supplementary materials). Math Biosci 283:13-29. https://doi.org/10.1016/j.mbs.2016.10.004

   [2] Arceo C, Jose E, Lao A, Mendoza E (2017) Reactant subspaces and kinetics of chemical reaction networks. J Math Chem 56(5):395-422. https://doi.org/10.1007/s10910-017-0809-x

   [3] Soranzo N, Altafini C (2009) ERNEST: a toolbox for chemical reaction network theory. Bioinform 25(21):2853-2854. https://doi.org/10.1093/bioinformatics/btp513