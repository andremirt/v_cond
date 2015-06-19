---
layout: page
title: Run
permalink: run/
---

The computation of \\(v_{cond}(\mathbf{r})\\) is done in two major steps

<ul>
  <li>compute a CI wavefunction GAMESS-US for your molecule</li>
  <li>use the CI wavefunction to compute \\(v_{cond}(\mathbf{r})\\) with vcond.exe</li>
</ul>

To demonstrate the usage of the v_cond package we give an example for the Hydrogen molecule. 

<h2>Run GAMESS-US</h2>

We compute a FullCI wavefunction in an aug-cc-pVTZ basis set. Use the GAMESS-US input file <code>h2_aug-cc-pVTZ_FCI.inp</code> from <a href="https://github.com/andremirt/v_cond/blob/master/examples/h2_aug-cc-pVTZ_FCI/">v_cond/examples/h2_aug-cc-pVTZ_FCI</a> or create it with the content

<pre>
<code  class="language-bash" data-lang="bash"> $CONTRL
   SCFTYP=RHF
   CITYP=ALDET
   RUNTYP=ENERGY
   ISPHER=1
!   EXETYP=CHECK
 $END
 $SYSTEM
   TIMLIM=1000000
   mwords=100
 $END
 $CIINP
   NRNFG(6)=1
 $END$
 $CIDET
   NCORE=0
   NACT=45
   NELS=2
 $END
 $BASIS
   GBASIS=ACCT
 $END
 $GUESS
   GUESS=HUCKEL
 $END
 $DATA
H2 RHF/aug-cc-pVTZ+FULLCI exp. geom
Dnh 4

HYDROGEN 1.0 0.0 0.0 0.373
 $END</code>
</pre>

Execute GAMESS-US
<pre>
<code  class="language-bash" data-lang="bash">rungmsvcond h2_aug-cc-pVTZ_FCI.inp > h2_aug-cc-pVTZ_FCI.inp</code>
</pre>

<h2>Run vcond.exe</h2>

First we process the GAMESS-US output to create the vcond.exe input by using the <code>cinput.sh</code> script. This script can be found in the v_cond source directory. If desired, link <code>cinput.sh</code> to your <code>~/bin</code> folder for easy execution
<pre>
<code  class="language-bash" data-lang="bash">ln -s ~/path_to_v_cond/v_cond/scripts/cinput.sh ~/bin/.</code>
</pre>

Consider to add execute permissions to <code>cinput.sh</code>
<pre>
<code  class="language-bash" data-lang="bash">chmod u+x ~/path_to_v_cond/v_cond/scripts/cinput.sh</code>
</pre>

Now we are ready to create the vcond.exe input by calling
<pre>
<code  class="language-bash" data-lang="bash">cinput.sh h2_aug-cc-pVTZ_FCI.out</code>
</pre>

This will create the files <code>vcond.basis</code> and <code>vcond.input</code>. 
Additionally you will need to provide the file <code>vcond.grid</code> with the points at which $$v_{cond}(\mathbf{r})$$ needs to be computed. For our example you can find a <code>vcond.grid</code> <a href="https://github.com/andremirt/v_cond/blob/master/examples/h2_aug-cc-pVTZ_FCI/">v_cond/examples/h2_aug-cc-pVTZ_FCI/</a>, to compute \\(v_{cond}(\mathbf{r})\\) along the bond axis.

In the final step we call
<pre>
<code  class="language-bash" data-lang="bash">vcond.exe</code>
</pre>
The files <code>vcond.dat</code>, <code>density.dat</code>, <code>rhovcond.dat, </code><code>vhartr.dat</code> and <code>wXC.dat</code>. <code>wXC.dat</code> contains the exchange-correlation energy density in the gauge of the exchange-correlation hole. See J. Chem. Theory Comput. 8, 3097 (2012) [<a href="../downloads/MirSeiGor-JCTC-12.pdf">pdf</a>] [<a href="http://dx.doi.org/10.1021/ct3003892">doi</a>] for more details.


