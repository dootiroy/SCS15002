Bibliography style (spr-chicago.bst) novelties.
Now, there is a possibility to choose a citation style: 
alphanumerical (e.g., [Ru01]), 
name-year (e.g., (Rudin, 2001) or Rudin (2001)) or number (e.g. [5]).
You can select it by writing an entry @SETTINGS to the bib file. That entry should 
have a field - options, with a value of citation style.
Default citation style for NALA journal is name-year (no options required).


Available options:
  alpha    -- alphanumerical citation style
  number   -- number citation style 
  nameyear -- name-year citation style
  unsort   -- turn off sorting (reference list is sorted by default) 
  nodoi    -- do not print doi (doi is printed by default)

Example 1: Choosing a number citation style.
  @settings{label,
     options = "number"
  }

Example 2: Choosing an alphanumerical citation style and not sorting references.
  @settings{label,
     options = "alpha,unsort"
  }

label can be any word, but all references labels have to be unique!


Important !!!
You need to write \nocite{settings label} 
in your tex file, otherwise  @settings will not become effective (unless you use
\nocite{*}, which will make a list of all the references appearing in bib file(s))

Example 3: Turning on @settings via tex file.
  bibfile.bib
  @settings{some_label, ....}

  texfile.tex
  ....
  \nocite{some_label}
  \bibliographystyle{spr-chicago}
  \bibliography{bibfile}



