function MedianFreq = Compute_Median_Frequency(norm2_s1, Wave_FreqS)
% norm2_s1 = abs(norm);
% Wave_FreqS = f;

MedianFreq = nan(size(norm2_s1,2),1) ;

for iTime = 1:size(norm2_s1,2)
    temp = norm2_s1(:,iTime);
    if sum(isnan(temp))==0 
        totinteg = trapz(temp(~isnan(temp))) ;
        Fmedsup = find(cumtrapz(temp(~isnan(temp)))>0.5*totinteg) ;
        MedianFreq(iTime,1) = Wave_FreqS(Fmedsup(1)) ;
    else
        toto = 1 ;
    end
end

