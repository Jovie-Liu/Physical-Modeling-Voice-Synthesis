# Instrument Modeling, Voice Synthesis, and Real-time Sound Effects

## Automatic FM Parameter Matching for Instruments

For any given instrument sounds, a reverse engineering algorithm is developed in MATLAB App Designer to search for the optimal FM parameters to resynthesize this instrument using FM formula. This method not only works for pitched instruments but also for percussions [[Presentation](https://drive.google.com/file/d/1dFm8WaWgDdPmAcj_LMZiInKogLoaoQlJ/view?usp=sharing)].

<img src="FM Parameter Matching/flute.png" style="width:800px">
<caption><center> Figure 1. Interface for FM parameter matching </center></caption>

## Digital Waveguide Physical Modeling of String Instruments

Here we use 1-dimensional [digital waveguide model](https://ccrma.stanford.edu/~jos/pasp/Digital_Waveguide_Models.html) to simulate the sounds of plucked string and struck string instruments. We tuned the model with recorded piano and pizzicato violin to demonstrate the model efficacy [[Video Demo](https://www.youtube.com/watch?v=zictMqwc3wc)].

<img src="Instrument Waveguide/waveguide.jpg" style="width:800px">
<caption><center> Figure 2. 1-D Waveguide Model </center></caption>

## Glottal Source - Vocal Tract Filter Modeling of Human Voice

Human voices are synthesized via a source-filter model where the derivative glottal source is generated from the LF-model and the vowel sounds are shaped and articulated via a piecewise-cylindrical vocal tract filter. The parameters in this model are highly tunable, which gives various qualities and flavors for the generated vowel sounds [[Video Demo](https://www.youtube.com/watch?v=PseuU-1j-qY)].

References:

[Glottal source modeling for singing voice synthesis](https://quod.lib.umich.edu/cgi/p/pod/dod-idx/glottal-source-modeling-for-singing-voice-synthesis.pdf?c=icmc;idno=bbp2372.2000.186;format=pdf)

[On the Transfer Function of the Piecewise Cylindrical Vocal Tract Model](http://musicweb.ucsd.edu/~trsmyth/pubs/smc21.pdf)

## Formant-Wave-Function Voice Synthesis and Real-time Sound Effects UI

The formant-wave function singing voice synthesis developed by IRCAM is modified into a real-time implementation with tunable parameters under which the qualities of the original voice could be modified to a great extent, such as the frequency range, vocal tone, pitch, relative amplitude of each spectral bands, dispersion and vibrato, etc. [[Video Demo](https://www.youtube.com/watch?v=tBkFWNNM0Uk)].

<img src="Singing Voice DAFX Real-time UI/UI.jpg" style="width:800px">
<caption><center> Figure 3. User Interface for Voice Modulation</center></caption>

Reference: [Time-Domain Formant-Wave-Function Synthesis](https://www.jstor.org/stable/3679809?seq=1)
