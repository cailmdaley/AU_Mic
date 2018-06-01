% AU Mic Collaborator Comments
% May 31, 2018

# General Questions & Comments

## MW

- Thanks again for including me on this. Having read the paper I see that including the results of my method to get a radial profile independent of assumptions about profile would involve quite a bit of change, particularly if the method needs to be described. There is also more work needed on our side to get a publication quality profile, most notably because my student's project was not concerned with details of the data like removing the flare, which could affect our conclusions. I suggest that your paper instead provides a great starting point for analysis with this method, and it would be great to have your carefully reduced data for such an analysis, and to include you and Cail on such a paper in due course.

### Typos: 

- [ ] 1 para 1 “knowN”, 
- [ ] 1 para 4 “Moshir & et al.|”, 
- [ ] page 6 note that Fig. 1 says “1.4mm” and section 3 “1.3mm”, 
- [ ] page 10 “asses”, 
- [ ] page 22 “aN Keplerian”


# Title

## MP

- [ ] The word "Hidden" doesn't seem like the best choice to convey that this mass is responsible for exciting the vertical structure. Not sure what to suggest-- perhaps something more clinical like "The Mass of Stirring Bodies in the AU Mic Debris Disk inferred from Resolved Vertical Structure".






# Abstract



## DW

- [ ] "smaller particles" suggests smaller than something; maybe "small particles"?
- [ ] "large planetesimals", Hilke already made this point, but I think this should express the range from a collection of large bodies, up to an Earth-size planet.



## MW

- [ ] I was concerned about the conclusion that you can rule out a gas giant or Neptune analogue, and when I got to 5.3 I think this is not right (more below).



## HS

- [ ] In the second line of the abstract, I suggest replacing ‘planetesimal bodies’ by ‘planetary bodies’.
- [ ] In the last line of the abstract, I suggest adding after ‘large planetesimals’ ‘or an earth-sized planet’






# 1. Introduction 



## HS

- [ ] In line 7: know - > known



## DW

- [ ] suggest to drop "(Pluto-sized)" here; one of the later the arguments (S5.3) suggests sizes might not be more than a few 100's km.
- [ ] for the reference to Hughes et al. 2018, suggest "see the review by Hughes et al. 2018 and references therein", to acknowledge better that there is a lot of literature on this topic.
- [ ] "the vertical structure can serve as a probe" sentence could use a citation, perhaps Thebault 2009.
- [ ] when using "inclination" and "eccentricity", these should be linked to the "orbits" of the dust, and not just "dust" (throughout the draft). Maybe the shorthand is fine since the meaning is not ambiguous, but better to be more precise with this language, I think.
- [ ] comma before "edge-on inclination"
- [ ] "4-6$\sigma$"? either pick one, e.g.  6 (as discussed in S4.3) and add something like "for our fiducial model", or use both but also indicate that this significance is model dependent?




# 2. Observations

## DW

- [ ] suggest to mention the wavelength before the longest baseline, since the wavelength should be known to discuss angular scales
- [ ] Table 2 would be easier to parse the flare if the exponents on the values were all the same, e.g. 10^2 microJy, or mJy. Maybe include a plot of the flare (flux vs. time), too?
- [ ] Regarding the flare, did you look at the polarization? In the Prox Cen data, there was a clear signature of polarization in a difference in XX and YY vs. time. It's not worth getting distracted by the flare, and maybe this analysis doesn't belong at all in this paper, or maybe it could be in an appendix if you have done the fits for XX and YY already.
- [ ] instead of "flare fits", how about "point-source fits"?
- [ ] "We exclude from our analysis [of the disk] the seven minutes..."
- [ ] Figure 1, consider to separate the panels in order to make it easier for readers to lift the left panel only from the paper, to show in talks


## SA

- [ ] In the observations, there is a typo claiming these are Band 7 data (top of p4).  There is also the claim that you have 2 GHz BW in the TDM windows: only 1.875 GHz is usable for those windows too (the 2 GHz is just a quirk of the OT).
- [ ] You claim to use MIRIAD (p4) in the calibration, but it doesn't seem like you do until you measure some emission.  My hope is it was only used for that!  Even so, why not just use the CASA viewer to do the same thing?
- [ ] The calibration discussion makes no mention of self-calibration.  Was that tried?  If not, why not?  If so, and it didn't help, it would be good to say so and report the interval(s) attempted.  [With a peak SNR ~ 20, I would think you could get a boost from self-cal.]


## MW

- [ ] Page 5: I think it is better to compare how the centroid shifted from when the flare was excluded to the epoch only including the flare, rather than to note how the centroid changed when excluding the flare (because you have already shown that it should be excluded).







# 3. Results
 
## MP

- [ ] Hilke mentioned the lack of significant brightness differences between the ansae and implications for the disk eccentricity (beginning of section 3). You have some info on the disk eccentricity from that but it depends on the disk's argument of periastron (if the argument of periastron is in the sky plane, the eccentricity would be the ratio of the fluxes minus unity, or 0.03 +/- 0.06, so consistent with zero). I could make a plot of the implied eccentricity as a function of the argument of periastron if you want to include a discussion of that - it would be something like "there is a 68% (or 95%, or 99%) chance the eccentricity is less than <a certain value>".


## SA

- [ ] When you quote an integrated line flux [p9], its a good idea to specify the *area* used to arrive at that number.


## DW

- [ ] "two limbs of the disk" doesn't sound right to me, maybe "two sides"?
- [ ] Figure 2, exactly how is "0" defined in the middle panel?
- [ ] suggest "features" instead of "clouds"
- [ ] recommend to give the gas limits its own subsection, to highlight this
result (and also make it easier to find).


## MW

- [ ] Fig. 2: What I am missing from this is a feeling for the typical uncertainties in these profiles, as well as how good a fit a Gaussian profile is to these vertical profiles.
- [ ] Page 9: I found it odd to see you talk of intensity maxima that are symmetric at ~10au, when Fig. 2 shows that the SE peak at that location is coincident with a minimum in the NW profile. I see you want this interpretation for later, but I think it needs to be described more clearly here. In fact, on page 15 you note that the maxima are not symmetric.



## HS 

- [ ] At the very top of page 7: I think there should be a comma after ‘Note’ - but I am not sure
- [ ] Also ion page 7, I suggest adding one sentence between ‘millimeter wavelength.’ and ‘Additional’ that explains what the absence of any significant brightness difference between the two limbs implies - and maybe even cite Margaret's paper on this. - Margaret, do you want to send one or two sentence about this?






# 4. Analysis

## MP 
- [ ]  could you mark the location of the star and/or the disk midplane in the residuals plots in Figure 5? I think that would make it more obvious which residual peaks are being discussed in the caption.


## DW

- [ ] A basic assumption that underlies the analysis is that a Gaussian profile in the vertical dimension describes the collection of dust orbits. This is not so obvious. Some justification from dynamics would be comforting. I know for our Kuiper Belt that the vertical distribution can be approximated very well by the sum of *two* Gaussians, see Brown 2001, AJ, 121, 2804.
- [ ] Table 3, if the Likelihood is going to reported here, then I think its definition should be in the text in S4.1.


## SA

- [ ] I hope you wouldn't get this, but I hear rumors of ApJ editors that are nitpicking on priors.  What you've done is fine (one could argue that choosing non-informative priors is not really Bayesian, but whatever), but there is a formally incorrect statement: you can't have a uniform prior with no bounds (its mathematically impossible), so just say you have some safely large bounds.
- [ ] On p11, I was confused by the statement about the inclination that "MCMC chain simply converges to two positions symmetric about 90 degrees."  What does that mean?
- [ ] How did you decide on the disk center?  Why not fit the center position too?
- [ ] I was a bit confused about the formulation of the model, particularly with respect to the vertical distribution.  You define h = H(r) / r, but maybe it would pay to explicitly state that h is independent of r.  Presuming that's correct, why would h be independent of r?  Maybe there's a theoretical reason, but its worth describing that.
- [ ] Speaking of which, is there a theoretical reason that the vertical distribution of particles would have a Gaussian density profile?  If so, a description of that would be helpful.


## MW

- [ ] Page 9: Can you call it “a CONSTANT aspect ratio”? As is I thought that H(r) is an arbitrary function of r whereas in fact it is linear because h is constant. Also, a small point is that you imply in 4 that the degeneracy leads you to adopt appropriate ray tracing methods. But here we are looking at optically thin thermal emission for which sophisticated ray tracing methods are not needed.
- [ ] Page 13: You say that h is uncorrelated with other aspects of the structure, but it is worth noting that h is correlated with inclination, because a disk that is edge-on and vertically broad has a similar vertical profile as one that is completely flat and inclined.
- [ ] Page 15: In 4.3 you are talking of the scale height at 40au and I wondered if it is important that Fig. 2 shows that the vertical height is independent of projected distance. This is probably not the case because at small projected separations we are seeing emission from much further out. I think you mention this later, but it might be worth commenting here on why you are only talking about 40au.



## HS

- [ ] On page 10, in the first line after equation 3: I recommend to define ‘PA’ here, unless it was done earlier in the manuscript and I missed it






# 5. Discussion



## SA

- [ ] It might be worthwhile to make a few flippant comments comparing the mm scale height for AU Mic to primordial disks.  For example, HL Tau seems to have an even lower aspect ratio (h ~ 0.01; see Pinte et al. 2016).  There's a few edge-on cases with similar, ~few percent, aspect ratios.  I certainly recognize the physical differences in these scenarios, but its interesting to think about this in an evolutionary context.



## MW

- [ ] Page 17: I couldn’t tell what the upper limit on gas mass refers to. Is this just the CO mass, or do you include some scaling for hydrogen assuming an abundance?
- [ ] S5.1:  
I think your conclusions are all valid. The inner edge of the disk is not so well constrained because there can be less emission here if the profile is steep, and part of the problem is that you have restricted the models to be power laws. This is fine for a first analysis, so I’m not suggesting to redo this, since it also motivates more complicated models such as the one I have been looking into that makes no assumption about the form of the profile. To forewarn you of the likely conclusions from that work: We did infer significant emission just inside 10au, but I don’t believe this yet because the method was affected by the stellar emission which moreover has the flare to contend with (but at least we don’t rule out a 10au ring like you are proposing). We got a slightly different profile on the two sides of the disk, but perhaps this difference will go when redoing the analysis with your carefully reduced data. In general though I would say that our derived profile agreed quite closely with the profile from Macgregor+2013, which is slightly different from yours. We also found the outer edge was not sharp, a possibility not explored in your analysis. I think these differences come down to your assumption of a single power law. One comment for this section is to pose the question of whether we can use the SED to distinguish between the bump at 10au projection being caused by emission at 10au or 40au? That is, would there be too much mid-IR emission from a ring at 10au at that level?
- [ ] Page 22: There is a discussion in section 4.9 of Marino et al. (2016) on using ALMA data to constrain the scale height of disks that are inclined to our line of sight. In fact that paper only gets an upper limit for the scale height in HD181327, but the discussion is relevant I think to this paper and this part of the discussion, since it was trying to motivate that this would be possible, as later shown in Kennedy et al. (2018) for HR4796.
- [ ] Page 23: I think it is incorrect to say that the collision velocities never exceed the escape velocities, since for some objects this can be the case. Perhaps you mean that “The mean collisional velocities vrel…”.
- [ ] Page 24: The Pan & Schlichting model assumes that the big bodies are embedded in the disk, i.e., they are the biggest bodies in the collisional cascade. There is no reason why there could not be a massive planet outside the disk which is causing stirring at this level. The planet does not need to be close to the disk to do this, since it can stir the disk through secular perturbations (e.g., Mustill & Wyatt 2009 consider eccentricity stirring, but a similar argument applies to an inclined planet). Of course the planet could be close to the disk, and in this case there probably would be constraints on how massive the planet is, though I’m less sure that the way you derive this is correct.
- [ ] Page 25: I agree with eq. 6, except that the “3" should be “(3+e)”, noting that this applies to the Hill radius at pericentre, and a different expression holds for the Hill radius at the planet’s apocentre. See Appendix B of Pearce & Wyatt (2014) for the derivation. Your concern about what happens when e goes to 1 is not founded, because a.(1-e) is the pericentre distance, which is what you know (presumably) rather than the semimajor axis. However, as mentioned above, I’m less convinced by how you are using this.



## HS

- [ ] Figure 6, in the figure legend: I suggest changing ‘Macgregor et al (2013)’ to ‘MacGregor et al. (2013)’.
- [ ] Also in Figure 6, I found the r < 15 au label confusing given the information in the figure caption and ‘r’ in the x-axis, do you mean ‘r_min < 15 au’?
- [ ] In section 5.3: I suggest to reword the first sentence in the paragraph after equation 5 to something like: ‘The velocity dispersion of the dust grains will be excited to about the escape velocity of the largest bodies governing the disk dynamics in the absence of any significant damping of their velocity dispersion. This result arises because viscous stirring has a larger cross section (i.e. shorter characteristic timescale) than collisions as long as the velocity dispersion is less than the escape velocity of the largest bodies that dominate the stirring. The two cross sections (and timescales) become comparable as v_rel approaches v_esc and collisions start to dominate limiting the growth of v_rel to about  v_esc (e.g. Schlichting 2014).
- [ ] On page 24, I recommend adding in the first line of the text ‘total’ before ‘dynamical mass'
- [ ] On page 24, in the first sentence of the new paragraph, I suggest replacing ‘collisional velocities’ to ‘velocity dispersion'
- [ ] On page 25, I would use the hill radius definition without the eccentricity here - which is just a (mp/3M_sun)^(1/3), then you can also leave out figure 7.
- [ ] Also, we should change the second line on page 25 from ‘stirring the disk’ to ‘stirring a dynamically cold disk’ to emphasize that the velocity dispersion is assumed to be sub-hill here.



## EIC 

- [ ] I agree with Hilke that the eccentricity dependence of the Hill sphere radius in equation (6) should be suppressed (because it is being used out of context), and that Figure 7 should be dropped accordingly.



## DW

### S5.0

- [ ] The "dust mass" here traces back to the flux density and assumed opacity, so maybe it is really the flux density measurements that should be compared?.
- [ ] I think it would be good break out the discussion of the gas into a separate paragraph, as above.

### S5.1

- [ ] Figure 6, since all the fits are power-laws, would it make more sense for this log-linear plot into a log-log plot?
- [ ] "Neither [value of the inner radius] parameter is well constrained"
- [ ] I am confused by the paragraph that compares the scattered light to the millimeter.  The sentence "The micron-sized grains... confined to a narrow birth-ring at 40 au" doesn't sound right-- the narrow birth ring is the source of micron-sized grains that blow out into a halo of bound and unbound orbits. And I don't see how the last sentence follows. Maybe I am just missing something obvious here.

### S5.3

- [ ] The stellar mass is taken as 0.5 solar masses for the Keplerian velocity. Should an uncertainty on that value be considered in subsequent calculations?
- [ ] "The collisional velocities..." paragraph contains *two* key results; I think this could be more clear if the two arguments were made in separate paragraphs.









# 6. Conclusions

## DW

- [ ] 4$\sigma$, as above, be consistent about 4$\sigma$ or 6$\sigma$ (or explain)
- [ ] The discussion of a planetary origin for the possible ring at 10 AU would fit better in S5.1 on the radial structure, rather than in the Conclusions.
- [ ] "ALMA Band 9", also give the wavelength (450 microns)






# Acknowledgements

- [ ] Finally, could you please add to the acknowledgments "MP gratefully acknowledges support from NASA grants NNX15AM35G and NNX15AK23G." ?
- [ ] Please add: H.E.S gratefully acknowledges support from NASA grant NNX15AK23G.
