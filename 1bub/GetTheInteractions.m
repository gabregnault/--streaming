function [INT] = GetTheInteractions(M,F,A,PHI)
test = zeros(length(M),length(F));
for ii = 1:length(M)
    for jj = 1:length(F)
        if abs(A(ii,jj)) > 0
            test(ii,jj) = 1;
        end
    end
end

s = 1;
for ii = 1:length(F)
    for jj = 1:length(M)
        for kk = jj:length(M)
            SUM = test(jj,ii) + test(kk,ii);
            if SUM == 2
                if M(jj) + M(kk) > 0
                    INT(s).order1 = M(jj);
                    INT(s).order2 = M(kk);
                    INT(s).Amp1 = A(jj,ii);
                    INT(s).Amp2 = A(kk,ii);
                    INT(s).freq = F(ii);
                    INT(s).Phase_shift = PHI(jj,ii) - PHI(kk,ii);
                    s = s+1;
                end
            end
        end
    end
end
end

