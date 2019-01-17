#  Input

## Accepted file formats

RMS currently supports its own `.rms` YAML based file format and Chemkin (`.inp`) format files (if an RMG species_dictionary.txt file is available it can be taken along with the Chemkin file to provide molecular structure information for the Chemkin file
  species).  Chemkin files can be in any units and individual units (quantities specified in Unitful format:  `1.0u"m^3"`) are accepted with the RMS format files.  These can be loaded

```
phaseDict = readinput("../src/testing/mech.rms")
```

or

```
phaseDict = readinput("chem_annotated.inp";
              spcdict="species_dictionary.txt")
```

respectively.  Note that `spcdict` is an optional parameter.  

## Output of Input File Reading

The output of reading an input file in RMS returns a dictionary of phase dictionaries indexed by the `name` of each phase.
Each phase dictionary has an array of Species objects corresponding to the key "Species" and an array of Reaction
Objects corresponding to the key "Reactions".

## YAML File Formats

RMS uses a YAML format input file.  YAML essentially stores information in embedded dictionaries and lists.
An example segment from a `.rms` file is given below.  
```
Phases:
- Species:
  - name: Ar
    smiles: '[Ar]'
    thermo:
      polys:
      - Tmax: 3459.6
        Tmin: 100.0
        coefs: [2.5, 9.24384602e-15, -1.36779837e-17, 6.66184769e-21, -1.00106912e-24,
          -1552.16105, 2.16745116]
        type: NASApolynomial
      - Tmax: 5000.0
        Tmin: 3459.6
        coefs: [2.49999999, 9.20455546e-12, -3.58608293e-15, 6.15198922e-19, -3.92041801e-23,
          -1552.16104, 2.16745122]
        type: NASApolynomial
      type: NASA
    type: Species
  - name: He
    smiles: '[He]'
    thermo:
      polys:
      - Tmax: 3459.6
        Tmin: 100.0
        coefs: [2.5, 9.24384602e-15, -1.36779837e-17, 6.66184769e-21, -1.00106912e-24,
          -1552.16105, -1.28349484]
        type: NASApolynomial
      - Tmax: 5000.0
        Tmin: 3459.6
        coefs: [2.49999999, 9.20455546e-12, -3.58608293e-15, 6.15198922e-19, -3.92041801e-23,
          -1552.16104, -1.28349478]
        type: NASApolynomial
      type: NASA
  name: phase
  Reactions:
  - kinetics: {A: 2.7590590000000007e-08, Ea: 26459.615999999998, n: 3.802, type: Arrhenius}
    products: [oxygen, octane]
    reactants: ['C[CH]CCCCCC', '[O]O']
    type: ElementaryReaction
  - kinetics: {A: 2.7590590000000007e-08, Ea: 26459.615999999998, n: 3.802, type: Arrhenius}
    products: [oxygen, octane]
    reactants: ['CCC[CH]CCCC', '[O]O']
    type: ElementaryReaction
```

Dashes `-` denote the beginning of a new `key:value` pair.  If done in series as above under `Species:` and under `Reactions:` this
makes an array of dictionaries.  Colons `:` denote the beginning of a `key:value` pair within a dictionary.  Key, value pairs following a dash that do not have a dash themselves are part of the same dictionary.  You can also define dictionaries and lists
as normal within Julia within the YAML.  

## .rms File Format

In YAML RMS amounts to a dictionary of phase dictionaries.  Each
phase dictionary has an entry corresponding to "Species", "Reactions" and "name".  The "name" corresponds to the name of the phase.  Within "Species" and "Reactions" are associated lists of dictionaries that correspond to Species and Reaction objects.  Beneath
"Species" and "Reactions" the structures all follow the same conventions.  

Each dictionary is assumed to correspond to an object within RMS that should be denoted by the `type` key and to have key:value maps
that correspond to all necessary parameters to construct the object corresponding to the `type` value.  Thus, the an Arrhenius
rate calculator can be defined as {A: 2.7590590000000007e-08, Ea: 26459.615999999998, n: 3.802, type: Arrhenius}.  In that case
RMS knows the object from the `type` value and what to put in the Arrhenius object from the remaining fields.  Of course for
more complex objects some of the values will correspond to objects themselves and thus be dictionaries with their own `type` values.  

Units can be defined for specific numerical values in the Unitful format `1.0u"m^3`.  However, currently molecular units are
not supported and thus moles must be used (this restriction applies only to `.rms` files and not to Chemkin format files).  
