% AU Mic Collaborator Comments on Version 2
% July 13, 2018

# Kevin 

- The title of section 3.1 should probably be 'CO content' instead of '18CO Content'?

- pg 12: 'The disk center is set by the measured position of the central point source, which were obtained in S2.' -> '..., which *was* obtained in S2.'

- caption of figure 6: 'A plot of the probability that the AU Mic disk has a given *eccentricity*'. In the second and third paragraph of section 5.3, do you mean the angle of periapsis (aka the orientation of the elliptical orbits) rather than the major axis? My understanding is that the major (and minor) axis are just based on the projected image of the disk, and hence by definition are in the plane of the sky. The angle of periapsis can vary from being in the plane of the sky to being along the line of sight, and having the angle of periapsis along the line of sight would produce a symmetric disk image regardless of the eccentricity. It also might help to add a sentence or two to the caption of figure 6 to explain how this probability is derived (e.g. 'The probability is based on the fraction of periapsis angles for a given eccentricity that are consistent with the observed SE-NW ansa flux density ratio (e.g. at low eccentricity any angle of periapsis would produce a small flux asymmetry, while at high eccentricity the angle of periapsis would need to be nearly aligned with the line of sight to produce small flux asymmetry)'.)


# Sean Andrews

- on p5, you discuss briefly the astrometric measurements being more precise than the measured stellar PM.  It would be good to quote some numbers here, since it sounds like a surprising statement given how large the stellar PM should be.  Even if you have more precision, are your positional measurements consistent with the GAIA-DR2 PM values?  If not, could the discrepancy be explained with phase noise in the ALMA data?

- The total flux density difference found wrt the MacGregor results (and extrapolations of previous single-dish measurements) is puzzlingly large.  You (rightly) point this out, but I never felt the paper offered a satisfactory explanation (or speculation) about *why* that's happened.  I'm especially confused in the context of Fig 7, since I would have assumed given the rough match in Sigma profiles that it would be hard to derive a factor of 2 different total flux density from such models.  What is going on here?

- Sect 3.1 title typo: (12CO, not 18CO).

- There's no perfect way to estimate a spectral line upper limit, but this is one case where you have a good handle on the geometry and a decent constraint on the host mass.  I think it would be useful to use a Keplerian mask constructed from this information + the spatial extent of the continuum, and use that to estimate the limit.

- Very minor point, but you specify a PA prior that is flat from 1-360 degrees.  The fact that you infer a PA posterior that has a single peak at PA = 128 degrees is formally inconsistent with that prior: it implies that your walkers didn't find what must be an equal-probability PA peak at 128+180 degrees (since you can't distinguish these from continuum data alone).  I would advocate specifying a flat prior from 1-180 degrees instead (since in practice I suspect this is what you've done).

- The idea of the "bumps" you've associated with the 10 au ring could alternatively be azimuthal asymmetries in the outer belt is compelling.  I suspect some dynamicists would be quite interested in an estimate of their radii and azimuthal angle separation in this context.

- I'm not exactly sure how it works, but its not super helpful to me to see a quoted significance (probability) on the AICc, but not on the delta(BIC).  How does one assess the meaning of a delta(BIC)?

- I would urge you to avoid ever saying that something is "beyond the scope of this work."  When I referee papers, I bristle at this (not infrequent) defensiveness.  Ideally, the referee has some input on what belongs or does not in a paper.  Its just polite to avoid shutting that down.  In this particular sentence (bottom of p24), I think this is unnecessary and could just be removed.  In other cases, the better approach would be to rephrase it...you know, like make a call for the community to consider doing that line of work, or whatever.

- Just to be clear, has Peter Plavchan given his consent for you to cite the AU Mic planet result as is done on p27?  If no, please ask!

- On p26, I'm confused about the "beam-subtracted quasar FWHMs" quoted...is that the deconvolved size of the calibrator emission?  If so, that seems very large to me!  At least in 1 dimension its quite similar to the inferred disk size.   Why are you so convinced "seeing" isn't a problem here?

- I got very confused on p28-29 about who "the author(s)" refers to (i.e., that its not self-referential).  I think this is just an attempt to avoid awkward over-citing...but it would be easier to read if you used the \citeauthor formatting designed for this kind of thing.

- The first sentence of Sect 6 (p32) brags about getting data at "double the angular resolution"...of course, you mean 2x *lower* angular resolution.  Doubling the resolution (i.e., the size of the synthesized beam) would be bad.


# Mark Wyatt 


- In 2 it would be worth saying which direction the offset is during the flare, as people will question if this is in the same direction as the fast moving features (they have done when I mentioned this to them and I didn’t have a reply).

- Fig. 2 plots the peak surface brightness. I found when modelling beta Pic in Telesco et al. 2005 and Dent et al. 2014 that it is better to use the flux integrated over the vertical extent, since that removes the possibility of the vertical size variations (including due to noise) and is a better reflection of how much emission there is as a function of distance. Hence this is what I use as the starting point for the modelling I was talking about to get the surface density distribution. (Actually, for beta Pic I used both the peak and the integrated flux, since then you can solve for both the surface density and the vertical height vs radius, and perhaps I can do the same with AU MIc). My recommendation would be to show both in Fig. 2, but it’s not the end of the world if you don’t.

- Fig. 5: Since you mention the second sub-plot being updated perhaps you are already addressing these issues, but I note that the colour scale does not make the 3sigma residuals clear as implied in the caption, and neither is it clear from c that the skinny model is a bad fit.

- Fig. 6 caption: I think you mean “has a given upper limit on its eccentricity”?

- In 5.3, you may find it helpful to refer to Fig. 7b of Wyatt et al. (1999), since that illustrates the point you are making about the degeneracy between the pericentre orientation and eccentricity to get a given asymmetry.

- * Since the abstract gives the lower limit as 1e-4Mearth, that should be in the text. For now it is usually quoted in radius, but also in kg and Mpluto. At the very end of 5.6, a similar calculation is performed, where the mass is given in Mearth, but there is a typo since it says 1e-3Mearth (should be 1e-4). There is also an issue of consistency since you are using 2g/cm^3 earlier, then 2.5g/cm^3 at the end. Not a big deal, but worth using the same I think. Actually, I see Eugene removed 1e-4Mearth from the abstract, so the comment is rather to make sure this is all consistent (and correct the typo).

- * I agree with Eugene that the discussion of the puffing up of the disk at shorter wavelengths could be tightened up. What I believe from this is that the observations in this paper give the vertical height of the parent belt needed for any modelling of this effect. Previous conclusions seem a bit all over the place, and will have made different assumptions. You do a good job of bringing these together, but since there is a range of short wavelength heights it’s hard to come to firm conclusions. What I take away from this is not that we find the mm height to be larger than the height at shorter wavelengths (contrary to expectations), since some estimates put this the other way round, rather that we don’t find strong evidence to support those expectations (which would have been the case say if the disk had been unresolved). We can also conclude that the vertical height in the short wavelength images is not set purely by radiation pressure as suggested could be the case in Thebault (2009). Eugene’s suggestion for why small grains may have a smaller vertical height is that these undergo higher damping due to higher collision rates, but I don’t think that is the answer, as they would undergo the same number of collisions before being destroyed and the higher rate just means this happens faster. On the other hand, damping could start to become more important if it is relatively harder to break up particles as they get smaller (which we know from the dispersal threshold that increases to smaller particles, but see also Krijt & Kama 2014). I think further work is needed to say what happens so I’m not convinced we should make strong statements here. Maybe wording like “These observations provide an important datapoint for the vertical height of large particles in the cascade from which to consider how this connects to the distribution of smaller particles. If the vertical height is found to decrease with particle size (as might be inferred from some estimates of the short wavelength vertical height) this could point to an increased role for damping with decreasing particle size” which could go on "perhaps due to the increased strength of such particles (e.g., Krijt & Kama 2014).”

- In the second last para of 5.6, typo “pertrubers”, and I think it should be “scale height IS EXPECTED TO DECREASE with decreasing radius”, since this is an assumption. I’m also not entirely convinced by this argument because it depends also on the surface density distribution. The vertical profile observed say in the direction of the star will be the weighted average of those at different radii. This would only be the size at the outer edge if most of the emission also comes from the outer edge.

- * The last para of 5.6 seems to be there to show that this does not contradict a model that may or may not be the correct explanation for the fast-moving features, which is more relevant for that paper rather than this one. I’m ok with that. But it would be helpful to have more information here, e.g. by adding something like “… in which the fast-moving features are explained by the ongoing destruction of a Varuna-sized progenitor”. I also didn’t get the argument that followed. Is the calculation that the lifetime is 3e4 yr and so with only one event we have a 3e4/2e7 chance of seeing this, which is increased if there are more events? Or does this account for the probability of these events occurring, since having 1000 Varuna’s still doesn't guarantee that we’d expect to see this. Perhaps this is all clear by a more thorough re-read of Chiang & Fung (2017), but that shouldn’t be necessary to understand the argument in this paper.

- * One quick calculation which may be worth including from the Mtot<1.5Mearth upper limit is to relate this to the Mcm=0.01Mearth in up to cm-sized grains derived from the mm-flux. Assuming the size distribution follows n(D) ~ D^-alpha up to Dmax, then Mtot / Mcm = (Dmax/Dcm)^(4-alpha), which means that for Dmax>420km requires alpha>3.71. This is not too far from that expected for particles in the strength regime (see eq. 24 of Wyatt, Clarke & Booth 2011), though this should only apply up to ~km in size, and the size distribution should become flatter above km in size as planetesimals become stronger due to self gravity. This could indicate that the size distribution is steeper than expected for a steady state collisional cascade if it is to extend up to the size of bodes required to stir the disk.
