      ******************************************************************
      ***** Range-dependent Acoustic Model, Version 1.5, 13-Sep-00 *****
      ******************************************************************

      This code was developed by Michael D. Collins at the Naval
      Research Laboratory in Washington, DC. It solves range-dependent
      ocean acoustics problems with the split-step Pade algorithm
      [M. D. Collins, J. Acoust. Soc. Am. 93, 1736-1742 (1993)]. A
      users guide and updates of the code are available via anonymous
      ftp from ram.nrl.navy.mil.

      Version 1.5 contains a correction to a bug in the dimension of
      quantities passed to subroutines fndrt and guerre that Laurie
      Fialkowski noticed.

      Version 1.4 contains a correction to a minor bug in subroutine
      guerre that Dave King noticed (amp1 and amp2 were declared
      twice) and a few other minor improvements.

      Version 1.3 contains a new root-finding subroutine.

      Version 1.2 contains a minor modification. The output to tl.grid
      is no longer zeroed out along the ocean bottom. This was done in
      previous versions so that the ocean bottom would be highlighted
      in graphical displays. The graphics codes ramclr, ramctr, and
      ramcc read in the bathymetry from ram.in and plot the ocean
      bottom directly.

      Version 1.1 contains two improvements:

      (1) An improved self starter. Stability is improved by using the
      factor (1-X)**2 instead of (1+X)**2 to smooth the delta function.
      The factor (1+X)**2 is nearly singular for some problems involving
      deep water and/or weak attenuation. Numerical problems associated
      with this singularity were detected by Eddie Scheer of Woods Hole
      Oceanographic Institute.

      (2) Elimination of underflow problems. A very small number is
      added to the solution in subroutine solve to prevent underflow,
      which can adversely affect run time on some computers. This
      improvement was suggested by Ed McDonald of the SACLANT Undersea
      Research Centre.

