Hi All,

 

Please see the attached screen shots.

 

Box plot show the distribution of gene classes.

 

The histogram shows the distribution of posterior probabilities for all genes. The essential (orange) and non-essential (pink) classes were fitted with a mixture of beta distributions

 

75% quantiles derived from the gold list (high essential and high non-essential) are represented by dashed black vertical lines (~0.21, ~0.97)

 

The scatter plot on top of the distribution shows the actual gold lists .

 

We could derive the cutoffs based on the quantiles of the two fitted beta distribution, Using same 75% for both give 0.2 and 0.99

 

Using different quantiles for each beta will make the cutoffs fall in visually pleasing places, but I cannot quite justify it.

 

Kourosh

From: Elsworth, Brendan <Brendan.Elsworth@fda.hhs.gov>
Date: Saturday, March 16, 2024 at 5:31 PM
To: Kourosh Zarringhalam <Kourosh.Zarringhalam@umb.edu>, Sida Ye <Sida.Ye001@umb.edu>, brendanelsworth <brendanelsworth@gmail.com>
Cc: Duraisingh, Manoj T. <mduraisi@hsph.harvard.edu>, Dass, Sheena <sdass@hsph.harvard.edu>
Subject: RE: [EXTERNAL] Transposon Manuscript thread -- March 10

Hey Kourosh,

 

Can you explain where the 75% comes from? And is it related to the fit or to the gold plus list?

 

To me, I would have thought the cut offs should include 95+% of the gold plus lists (maybe removing outliers – not sure about that). Seems like a pretty significant chunk – especially around the non-essentials are being missed by the cut off?

 

Thanks,

Brendan

From: Duraisingh, Manoj T. <mduraisi@hsph.harvard.edu>
Sent: Saturday, March 16, 2024 5:55:01 PM
To: Elsworth, Brendan <brendan.elsworth@fda.hhs.gov>; Kourosh Zarringhalam <Kourosh.Zarringhalam@umb.edu>; Sida Ye <Sida.Ye001@umb.edu>; brendanelsworth <brendanelsworth@gmail.com>
Cc: Dass, Sheena <sdass@hsph.harvard.edu>
Subject: Re: [EXTERNAL] Transposon Manuscript thread -- March 10

 

CAUTION: EXTERNAL SENDER

Looks very cool.

As an endpoint I think the kind of cut-offs/calls that would be helpful would result in groups for essential and non-essential genes that include about 30-35% of the total numbers of genes each, with very high confidence? Thoughts?

From: Kourosh Zarringhalam <Kourosh.Zarringhalam@umb.edu>
Sent: Saturday, March 16, 2024 6:38 PM
To: Duraisingh, Manoj T. <mduraisi@hsph.harvard.edu>; Elsworth, Brendan <Brendan.Elsworth@fda.hhs.gov>; Sida Ye <Sida.Ye001@umb.edu>; brendanelsworth <brendanelsworth@gmail.com>
Cc: Dass, Sheena <sdass@hsph.harvard.edu>
Subject: Re: [EXTERNAL] Transposon Manuscript thread -- March 10

 

CAUTION: This email originated from outside of the organization. Do not click links or open attachments unless you recognize the sender and know the content is safe.

 

The cut offs are based on box plots (gold list). If you look at high confidence ones only, the 75% were calculated based on those.

 

There is quite a bit of heterogeneity in high confidence non essentials. Removing outliers will change these numbers of course, maybe push them closer to what you expect d.

 

I think one way to proceed is to select quantiles based on the fitted beta distributions, and pick them so they would include a given proportion (either from gold list or all genes.)

 

Kourosh 

 

Get Outlook for iOS

From: Elsworth, Brendan <Brendan.Elsworth@fda.hhs.gov>
Date: Saturday, March 16, 2024 at 6:45 PM
To: Kourosh Zarringhalam <Kourosh.Zarringhalam@umb.edu>, Duraisingh, Manoj T. <mduraisi@hsph.harvard.edu>, Sida Ye <Sida.Ye001@umb.edu>, brendanelsworth <brendanelsworth@gmail.com>
Cc: Dass, Sheena <sdass@hsph.harvard.edu>
Subject: RE: [EXTERNAL] Transposon Manuscript thread -- March 10

CAUTION: EXTERNAL SENDER

I think the high confidence is too small a set. Even though I labelled them high and medium, they are both quite strict and I think they should be combined (they do not overlap).

From: Kourosh Zarringhalam <Kourosh.Zarringhalam@umb.edu>
Sent: Saturday, March 16, 2024 8:25 PM
To: Elsworth, Brendan <Brendan.Elsworth@fda.hhs.gov>; Duraisingh, Manoj T. <mduraisi@hsph.harvard.edu>; Sida Ye <Sida.Ye001@umb.edu>; brendanelsworth <brendanelsworth@gmail.com>
Cc: Dass, Sheena <sdass@hsph.harvard.edu>
Subject: Re: [EXTERNAL] Transposon Manuscript thread -- March 10

 

CAUTION: This email originated from outside of the organization. Do not click links or open attachments unless you recognize the sender and know the content is safe.

 

I combined the high and mid confidence gold list to produce two classes only. Using a bootstrap sampling (100 bootstrap), we can calculate the mean and sd of HMS for each class and use mean +/- 2 * sd as confidence interval cutoffs, which gives ~0.21 and ~0.82 cutoffs (attached).

 

These cutoffs result in ~35% essential and 43% non-essential genes overall.

 

This also results in ~76% and 83% of gold list essential and non-essential genes to be correctly classified respectively.

 

Kourosh

From: Elsworth, Brendan <Brendan.Elsworth@fda.hhs.gov>
Date: Saturday, March 16, 2024 at 8:45 PM
To: Kourosh Zarringhalam <Kourosh.Zarringhalam@umb.edu>, Duraisingh, Manoj T. <mduraisi@hsph.harvard.edu>, Sida Ye <Sida.Ye001@umb.edu>, brendanelsworth <brendanelsworth@gmail.com>
Cc: Dass, Sheena <sdass@hsph.harvard.edu>
Subject: RE: [EXTERNAL] Transposon Manuscript thread -- March 10

CAUTION: EXTERNAL SENDER

It certainly all sounds and looks reasonable to me. From a statistical point of view, are you happy with this method or does it feel like we are stretching it to be what we want?

 

Brendan

From: Kourosh Zarringhalam <Kourosh.Zarringhalam@umb.edu>
Date: Saturday, March 16, 2024 at 10:18 PM
To: Elsworth, Brendan <brendan.elsworth@fda.hhs.gov>, Duraisingh, Manoj T. <mduraisi@hsph.harvard.edu>, Sida Ye <Sida.Ye001@umb.edu>, brendanelsworth <brendanelsworth@gmail.com>
Cc: Dass, Sheena <sdass@hsph.harvard.edu>
Subject: Re: [EXTERNAL] Transposon Manuscript thread -- March 10

WARNING: Harvard cannot validate this message was sent from an authorized system. Please be careful when opening attachments, clicking links, or following instructions. For more information, visit the HUIT IT Portal and search for SPF.

It makes sense from a statistical point of view. The quantile version for each beta is also statistically sound, but does not take into account training data (gold list) into picking quantiles. So it’s justified.

 

Kourosh

From: Duraisingh, Manoj T. <mduraisi@hsph.harvard.edu>
Sent: Sunday, March 17, 2024 9:17 AM
To: Kourosh Zarringhalam <Kourosh.Zarringhalam@umb.edu>; Elsworth, Brendan <Brendan.Elsworth@fda.hhs.gov>; Sida Ye <Sida.Ye001@umb.edu>; brendanelsworth <brendanelsworth@gmail.com>
Cc: Dass, Sheena <sdass@hsph.harvard.edu>
Subject: Re: [EXTERNAL] Transposon Manuscript thread -- March 10

 

CAUTION: This email originated from outside of the organization. Do not click links or open attachments unless you recognize the sender and know the content is safe.

 

All looks very promising, Kourosh.

 

Brendan, could you send us your write-up of how the gold list candidates were chosen, listing the criteria?

 

For one, these could be limited to genes with >6 TTAAs for instance if the Pf transposon mutagenesis data was used.

 

Best,

Manoj

From: Elsworth, Brendan <Brendan.Elsworth@fda.hhs.gov>
Date: Sunday, March 17, 2024 at 10:40 AM
To: Duraisingh, Manoj T. <mduraisi@hsph.harvard.edu>, Kourosh Zarringhalam <Kourosh.Zarringhalam@umb.edu>, Sida Ye <Sida.Ye001@umb.edu>, brendanelsworth <brendanelsworth@gmail.com>
Cc: Dass, Sheena <sdass@hsph.harvard.edu>
Subject: RE: [EXTERNAL] Transposon Manuscript thread -- March 10

All filters are in the powerpoint. Full set of data is in the attached table. It has number of sites and insertions in Pf data.                                                                                                                                                                                                                                            

 

I thought about limiting it to number of sites but given that it either has to be consistent across all 3 screens (pf/Pb/Tg) or across Pf/Pb screens and both have targeted KOs, it is very stringent. Adding site/gene length would bias our data and not accurately show what fewer TTAA essentials looks like.

 

Best,

Brendan

From: Duraisingh, Manoj T. <mduraisi@hsph.harvard.edu>
Sent: Sunday, March 17, 2024 10:52 AM
To: Elsworth, Brendan <Brendan.Elsworth@fda.hhs.gov>; Kourosh Zarringhalam <Kourosh.Zarringhalam@umb.edu>; Sida Ye <Sida.Ye001@umb.edu>; brendanelsworth <brendanelsworth@gmail.com>
Cc: Dass, Sheena <sdass@hsph.harvard.edu>
Subject: Re: [EXTERNAL] Transposon Manuscript thread -- March 10

 

CAUTION: This email originated from outside of the organization. Do not click links or open attachments unless you recognize the sender and know the content is safe.

 

Not sure I agree with your point of the data being biased if you limit it in any way. I see that as stringency. I think we are trying to produce a list of essential proteins as a gold list with as much prior confidence as possible, before using it with the Pk data.

How many genes do you lose if you only include Pf genes with signficant numbers of TTAAs?

Best,

Manoj