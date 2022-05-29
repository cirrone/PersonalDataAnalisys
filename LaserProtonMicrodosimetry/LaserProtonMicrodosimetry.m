clear all
close all
Scorer = ReadPHSPFile_2('Scorer_Pablo.phsp');

% Particle arrival time
%
ParticleArrivalTime = Scorer.VarName10;

% Particle kinetic energy at the detector surface
%
ParticleKineticEnergy = Scorer.VarName6;

%[countsE_kin, binE_kin] = hist(Scorer.VarName6, 300);


%[count_time, bin_time] = hist(Scorer.VarName10, 300);

plot(ArrivalTime, ParticleKineticEnergy);

