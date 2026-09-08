function SpecEnt = Compute_Spectral_Entropy(norm2_s1, Wave_FreqS);

% norm2_s1=norm;

SpecEnt = nan(size(norm2_s1,2),1) ;

for iTime = 1:size(norm2_s1,2)
    SpecEnt(iTime,1) = Spectral_Entropy(norm2_s1(:,iTime),Wave_FreqS);
end
