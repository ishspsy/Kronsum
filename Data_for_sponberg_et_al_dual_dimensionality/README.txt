These data files are associated with the manuscript “Dual dimensionality reduction reveals independent encoding of motor features in a muscle synergy for insect flight control” by Simon Sponberg, Thomas L. Daniel, and Adrienne L. Fairhall. 

In general the data are taken from EMG recordings of the downstroke flight muscles (DLMs) of hawkmoths, Manduca sexta, flying in a visual arena while tethered to a torque-meter that measured turning torque in the yaw direction.

These data were originally collected for the following publication:
Sponberg, S., & Daniel, T. L. (2012). Abdicating power for control: a precision timing strategy to modulate function of flight power muscles. Proceedings of the Royal Society B-Biological Sciences, 279(1744), 3958–3966. doi:10.1098/rspb.2012.1085

The data are organized into many individual wingstrokes taken from 7 individual moths identified as “J”, “K”, “L”, “M”, “N”, “P”, and “Q”
Wingstrokes were delineated using either spike-triggering or phase-triggering methods as described in the manuscript. These produce ensembles of wingstrokes consisting of 500 samples of instantaneous torque measurements immediately following the triggering event (either the right DLM’s spike or the zero crossing of the phase variable). Since torque was samples at 10 kHz, this corresponds to 50 ms or the typical wingstroke duration during tethered flight. Each wingstroke also has associated EMG data consisting of spike times and other measures specified below.

The data files include five types:

Filetype #1: Tspiketrigensemble*.csv, where * is one of the individual moth labels (e.g. “J”). This is an N x 500 csv file where N is the number of wingstrokes and each of the 500 columns is one instantaneous torque measurement. These is the spike-triggered ensemble data.

Filetype #2: Tphaseensemble*.csv, where * is one of the individual moth labels (e.g. “J”). This is an N x 500 csv file where N is the number of wingstrokes and each of the 500 columns is one instantaneous torque measurement. These is the phase-triggered ensemble data.

Filetype #3: wingstroke_data_*.csv, where * is one of the individual moth labels (e.g. “J”). This is an N x 15 csv file where N is the number of wingstrokes and each of the 15 columns are a specific measurement for each wingstroke. Each row corresponds to the same wingstroke in the ensemble files. Columns are as follows:

Note that time units are in seconds unless otherwise specified, the phase variable ranges from -pi to pi radians, missing data are ‘NaN’. 

1:    'animalnum' — a numeric animal identifier
2:    'trialnum' — which flight bout the wingstroke is from
3:    'Rightspiketime' — absolute time in the flight bout of the spike in the right DLM during that wingstroke
4:    'deltat' — timing difference between the left and right DLM spikes
5:    'right ISI' — the interspike interval between the right DLM spike in previous and current wingstroke
6:    'left ISI' — the interspike interval between the left DLM spike in previous and current wingstroke
7:    'right phase' — the phase at which the right DLM spike occurs.
8:    'left phase' — the phase at which the left DLM spike occurs.
9:    'stimpos' — normalized position (-1 to 1) of the sinusoidal visual stimulus at the beginning of the wingstroke (i.e. at the triggering event).
10:    'stimvel' — normalized velocity (-1 to 1) of the sinusoidal visual stimulus at the beginning of the wingstroke (i.e. at the triggering event).
11:    'original tor imp (V * ms)’ — Uncalibrated torque meter impulse in (V * ms) over the wingstroke. Mean torque-meter reading (in V) during the wingstroke is this value divided by 50 ms.
12:   'New tor imp (mN * m * s)’ — Calibrated torque impulse in (mN * m * s) over the wingstroke calibrated as described in the manuscript. Mean torque (in mN*m) during the wingstroke is this value divided by 0.05 s.
13:   'phase period' — The time between zero phase crossing events from the previous wingstroke to the current one.
14:    'right spike time from zero phase' — The time between the zero phase crossing event and the time of the right DLM spike (t_R in the manuscript)
15:    'left spike time from zero phase' — The time between the zero phase crossing event and the time of the left DLM spike (t_L in the manuscript)

Filetype #4: Features_spike_tiggered.mat contains one Matlab structure for each moth (e.g. “J”) that includes the following fields:

1) STA — the spike triggered average torque waveform (500 samples)
2) STAscores — one score for each of N wingstrokes representing the score of the STA for that wingstroke (i.e. the projection of that wingstroke onto the STA).
3-5) Features structures for the SVD (single value decomposition), SVDcentered (single value decomposition of the ensemble with the STA removed), and the PLS (partial least squares) methods for feature analysis as described in the paper. The SVD an SVDcentered structures contain a “diagS” field, which is the SVD sigma diagonal matrix or “S” value for each of the first 20 dimensions. The “scores” field is an N by 20 matrix giving the score for each of N wingstrokes for the first 20 features (dimensions). The ‘features’ matrix is a 500x20 matrix that defines the shape of the feature over the 500 sample long wingstroke for each of the first 20 features. In the PLSfeatures structure, “X” refers to the torque matrix (matrix M in the manuscript) and “Y” refers to the motor commands (matrix U in the manuscript). “XLs” and “YLs” define the X and Y loadings which are the shape of the feature for the first 20 dimensions. “XSs” and “YSs” define the scores for these features for each of N wingstrokes. “PCTVARs” is a 2x20 matrix, where the 2 rows are t_R and t_L respectively and each column is the percent variation in t_R or t_L that is explained by that feature. 

When data necessary for feature analysis were missing, the corresponding wingstroke was dropped from the Features data. Therefore the number of wingstrokes may not be exactly equal to the number present in the ensemble or wingstroke datasets above (e.g. Animal “J” had 547 wingstrokes, but only 535 were complete).

Filetype #5: Features_phase_triggered.mat is exactly like Features_spike_triggered.mat, but used the phase-triggered ensembles for the feature analysis rather than the spike triggered ensembles. STA fields in this case refer to the phase-triggered average.
