function Markov_chain_DNA_replication

% This function seeks to demonstrate the evolutionary advantage of left-right
% symmetry breaking in the kinetic parameters of DNA replica construction
% on top of a template strand. The competition between the requirement of
% low kinetic barrier for monomer induction and the high kinetic barrier
% for monomer retention leads to the symmetry breaking. 

% There are two time-scales in the problem: The rate of
% H-bonding/dissociation between the free monomers and the template strand,
% and the rate of covalent bond formation between monomers on the replica
% strand. When the monomer supply is abundant and consequently, the
% H-bonding rates are high, there is no evolutionary pressure. But when the
% supply is scarce, the H-bonding rates are low and the ability to retain
% attached monomers becomes the deciding factor for successful replication.
% The DNA replication is modeled as a Markov chain below. '0' represents
% the absence of a monomer and '1', its presence. The script models the
% growth of a 5-nt long template. Covalent bond formation is not explicitly
% included in the chain. It is factored in indirectly as the
% retention-advantage factor below, calculated through the probability for
% the chain to stay in the '11111' state. 

% The DNA replica strand construction is assumed to take place
% cooperatively, with the neighboring monomers hydrogen-bonded to the
% template influencing the rate of monomer attachment. This rate is not
% symmetric w.r.t the left and right neighbors, and can be different. This
% difference is encoded in the variable 'a' below, which can be varied from
% zero (infinitely high kinetic barrier) to infinity (instantaneous
% formation of H-bond). The output produced by the script would show that
% it is evolutionarily advantageous for the DNA heteropolymer to have its
% left-right symmetry broken, with low barriers to the right (left) for
% rapid induction of monomers and high barriers to the left (right) for
% retention of the monomers attached to the template. The variable names
% are nearly self-explanatory. 


% For more information on the model, please refer to "H.  Subramanian, 
% R. A. Gatenby, 'Evolutionary advantage of directional symmetry breaking 
% in self-replicating polymers', Journal of Theoretical Biology, vol. 446,
% pp. 128–136 (2018)".

breakage_rate = 0.5; % per sec. Hydrogen bond breakage rate. 
formation_rate = 1; % hydrogen bond formation rate, Per second. % Rate constant is usually 10^7 M^-1 Sec^-1 
h_freebonding_rate = formation_rate; % Rate of bond formation between two free-floating monomers. 
p_formation_rate = 10; % per sec. Rate of covalent bond formation. 


numdiv=30; % For discretization of the cooperativity parameter. 
a = 1:(4-1)/(numdiv-1):4;  % Cooperativity parameter, to raise or lower the kinetic barriers.  
b=a;

xticklabel = 1:0.5:4;
yticklabel = xticklabel(2:end);

numctrs=30; % number of contour lines


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 5-bond polymer growth dynamics through Markov Chain model:


xi = a; % rates of right hyd-bonding for different barrier heights
yi = xi; % rates of left hyd-bonding for different barrier heights
r0 = p_formation_rate/formation_rate; % Relative rate of formation of covalent bonds w.r.t hydrogen bond formation rates. 
rg = h_freebonding_rate/formation_rate; % Rate of processes that rob monomers away from polymerization, such as dimerization.


for inda=1:length(a)  % For various left-right asymmetric cooperativity factors, calculate the rates
    for indb=1:length(b)

            p01 = 1; q01 = breakage_rate/formation_rate; 
            pr1 = formation_rate*a(inda)/formation_rate; 
            qr1 = breakage_rate*a(inda)/formation_rate;
            pc1 = formation_rate*(a(inda)*b(indb))/formation_rate; 
            qc1 = breakage_rate*(a(inda)*b(indb))/formation_rate;
            pl1 = formation_rate*b(indb)/formation_rate; 
            ql1 = breakage_rate*b(indb)/formation_rate;

%           This function returns back the Markov chain transition matrix. 
%           If circular strand matrix is needed, use _circ. For linear strands, use _lin.             
            Aasym5 = generator5_circ(p01, q01, pr1, qr1, pl1, ql1, pc1, qc1); % Change circ and lin below as well.

            % Residence time at 11111 is the reciprocal of the diagonal element at
            % 11111
            residencetime5_asym(inda,indb) = abs(1/Aasym5(32,32));

            % hitting time for double-strand formation
            % Least square non-negative solution
            [htasym5(:,inda,indb), resnorm2] = lsqnonneg(Aasym5(1:31,1:31), -ones(31,1));

            if (resnorm2 > 10^-10) 
                warning('LSQ not working');
            end
            
            % This gives the probability that the chain stays in the
            % '11111' state
            pr5_C1(inda,indb) = 1-exp(-r0*residencetime5_asym(inda,indb)); 

            % This is the advantage of bonding with the template over
            % bonding with other free monomers.
            pr5_G1(inda,indb) = (1/htasym5(1,inda,indb))/(rg + (1/htasym5(1,inda,indb)));
            
%             residencetime5(inda,indb) = residencetime5_asym(inda,indb);
%             zippingtime5(inda,indb) = htasym5(1,inda,indb);

            
    end
end

% For comparison, evaluate the above factors for the non-cooperative case
% as well. The plots below exhibit the ratio of
% cooperative-over-noncooperative rate advantages. 
Anocoop5 = generator5_circ(p01, q01, p01, q01, p01, q01, p01, q01); % Change circ and lin here as well. 
residencetime5_nocoop = abs(1/Anocoop5(32,32));
[htnocoop5, resnorm1] = lsqnonneg(Anocoop5(1:31,1:31), -ones(31,1)); 
pr5_C3 = 1-exp(-r0*residencetime5_nocoop);
pr5_G3 = (1/htnocoop5(1))/(rg + (1/htnocoop5(1)));



figure; 

[cs,h]=contourf(xi,yi,(pr5_G1.*pr5_C1)./(pr5_G3.*pr5_C3),numctrs);
% clabel(cs,h,'LabelSpacing',150);
ax=gca;
axis square
set(h,'LineStyle','none');
xlabel('Rate modulating factor \alpha_L','FontSize',15);
ylabel('Rate modulating factor \alpha_R','FontSize',15);
set(ax,'YTick',yticklabel);
set(ax,'XTick',xticklabel);
set(ax, 'FontSize',15);
%  title('Probability for growth and C-bonding-Asym');
[cmin1,cmax1]=caxis;
c=colorbar;
set(c,'FontSize',15);
xlabel(c,'P/P_0','FontSize',15);

figure;
[cs,h]=contourf(xi,yi,pr5_G1./pr5_G3,numctrs);
% clabel(cs,h,'LabelSpacing',150);
ax=gca;
axis square
set(h,'LineStyle','none');
xlabel('Rate modulating factor \alpha_L','FontSize',15);
ylabel('Rate modulating factor \alpha_R','FontSize',15);
set(ax,'YTick',yticklabel);
set(ax,'XTick',xticklabel);
set(ax, 'FontSize',15);
%  title('Probability for growth');
[cmin,cmax]=caxis;
colorbar;
c=colorbar;
set(c,'FontSize',15);
xlabel(c,'P_g/P_{g0}','FontSize',15);

figure;
[cs,h]=contourf(xi,yi,pr5_C1./pr5_C3,numctrs);
% clabel(cs,h,'LabelSpacing',150);
ax=gca;
axis square
set(h,'LineStyle','none');
xlabel('Rate modulating factor \alpha_L','FontSize',15);
ylabel('Rate modulating factor \alpha_R','FontSize',15);
set(ax,'YTick',yticklabel);
set(ax,'XTick',xticklabel);
set(ax, 'FontSize',15);
%  title('Probability for C-bonding-Asym');
[cmin,cmax]=caxis;
colorbar;
c=colorbar;
set(c,'FontSize',15);
xlabel(c,'P_c/P_{c0}','FontSize',15);

disp('5-bond Nocoop fitness is');
disp(pr5_G3.*pr5_C3);

end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function A = generator5_lin(p0, q0, pr, qr, pl, ql, pc, qc)


% The order of the elements are as follows:
% 00000(1)
% 10000(2) 01000(3) 00100(4) 00010(5) 00001(6)
% 11000(7) 01100(8) 00110(9) 00011(10) 10001(11)
% 10100(12) 01010(13) 00101(14) 10010(15) 01001(16)
% 11100(17) 01110(18) 00111(19) 10011(20) 11001(21)
% 11010(22) 01101(23) 10110(24) 01011(25) 10101(26)
% 11110(27) 01111(28) 10111(29) 11011(30) 11101(31)
% 11111(32)


% group 1
basis(1,1:5) = '00000';
% group 2 (2-6)
basis(2,1:5) = '10000'; basis(3,1:5) = '01000'; basis(4,1:5) = '00100'; basis(5,1:5) = '00010'; basis(6,1:5) = '00001';
% group 3 (7-11)
basis(7,1:5) = '11000'; basis(8,1:5) = '01100'; basis(9,1:5) = '00110'; basis(10,1:5) = '00011'; basis(11,1:5) = '10001';
% group 4 (12-16)
basis(12,1:5) = '10100'; basis(13,1:5) = '01010'; basis(14,1:5) = '00101'; basis(15,1:5) = '10010'; basis(16,1:5) = '01001';
% group 5 (17-21)
basis(17,1:5) = '11100'; basis(18,1:5) = '01110'; basis(19,1:5) = '00111'; basis(20,1:5) = '10011'; basis(21,1:5) = '11001';
% group 6 (22-26)
basis(22,1:5) = '11010'; basis(23,1:5) = '01101'; basis(24,1:5) = '10110'; basis(25,1:5) = '01011'; basis(26,1:5) = '10101';
% group 7 (27-31)
basis(27,1:5) = '11110'; basis(28,1:5) = '01111'; basis(29,1:5) = '10111'; basis(30,1:5) = '11011'; basis(31,1:5) = '11101';
% group 8
basis(32,1:5) = '11111';

A(1:32,1:32)=0; % NON-CIRCULAR STRAND
A(1,2:6)=p0;

A(2,1)=q0; A(2,7)=pr; A(2,11)=p0; A(2,12)=p0; A(2,15)=p0;
A(3,1)=q0; A(3,8)=pr; A(3,7)=pl; A(3,13)=p0; A(3,16)=p0;
A(4,1)=q0; A(4,9)=pr; A(4,8)=pl; A(4,14)=p0; A(4,12)=p0;
A(5,1)=q0; A(5,10)=pr; A(5,9)=pl; A(5,15)=p0; A(5,13)=p0;
A(6,1)=q0; A(6,11)=p0; A(6,10)=pl; A(6,16)=p0; A(6,14)=p0;

A(7,2)=qr; A(7,3)=ql; A(7,17)=pr; A(7,21)=p0; A(7,22)=p0;
A(8,3)=qr; A(8,4)=ql; A(8,18)=pr; A(8,17)=pl; A(8,23)=p0;
A(9,4)=qr; A(9,5)=ql; A(9,19)=pr; A(9,18)=pl; A(9,24)=p0;
A(10,5)=qr; A(10,6)=ql; A(10,20)=p0; A(10,19)=pl; A(10,25)=p0;
A(11,6)=q0; A(11,2)=q0; A(11,21)=pr; A(11,20)=pl; A(11,26)=p0;

A(12,2)=q0; A(12,4)=q0; A(12,17)=pc; A(12,24)=pr; A(12,26)=p0;
A(13,3)=q0; A(13,5)=q0; A(13,18)=pc; A(13,25)=pr; A(13,22)=pl;
A(14,4)=q0; A(14,6)=q0; A(14,19)=pc; A(14,26)=p0; A(14,23)=pl;
A(15,5)=q0; A(15,2)=q0; A(15,20)=pr; A(15,22)=pr; A(15,24)=pl;
A(16,6)=q0; A(16,3)=q0; A(16,21)=pl; A(16,23)=pr; A(16,25)=pl;


A(17,7)=qr; A(17,8)=ql; A(17,12)=qc; A(17,27)=pr; A(17,31)=p0;
A(18,8)=qr; A(18,9)=ql; A(18,13)=qc; A(18,28)=pr; A(18,27)=pl;
A(19,9)=qr; A(19,10)=ql; A(19,14)=qc; A(19,29)=p0; A(19,28)=pl;
A(20,10)=q0; A(20,11)=ql; A(20,15)=qr; A(20,30)=pr; A(20,29)=pl;
A(21,11)=qr; A(21,7)=q0; A(21,16)=ql; A(21,31)=pr; A(21,30)=pl;

A(22,7)=q0; A(22,15)=qr; A(22,13)=ql; A(22,27)=pc; A(22,30)=pr;
A(23,8)=q0; A(23,16)=qr; A(23,14)=ql; A(23,28)=pc; A(23,31)=pl;
A(24,9)=q0; A(24,12)=qr; A(24,15)=ql; A(24,29)=pr; A(24,27)=pc;
A(25,10)=q0; A(25,13)=qr; A(25,16)=ql; A(25,30)=pl; A(25,28)=pc;
A(26,11)=q0; A(26,14)=q0; A(26,12)=q0; A(26,31)=pc; A(26,29)=pc;


A(27,24)=qc; A(27,22)=qc; A(27,18)=ql; A(27,17)=qr; A(27,32)=pr;
A(28,25)=qc; A(28,23)=qc; A(28,19)=ql; A(28,18)=qr; A(28,32)=pl;
A(29,26)=qc; A(29,24)=qr; A(29,20)=ql; A(29,19)=q0; A(29,32)=pc;
A(30,22)=qr; A(30,25)=ql; A(30,21)=ql; A(30,20)=qr; A(30,32)=pc;
A(31,23)=ql; A(31,26)=qc; A(31,17)=q0; A(31,21)=qr; A(31,32)=pc;

A(32,27)=qr; A(32,28)=ql; A(32,29)=qc; A(32,30)=qc; A(32,31)=qc; 

A = A - diag(sum(A,2)); % Make the rows sum to zero. Conservation of probability flow.

for i=1:32
    for j=i:32
        if sum(abs(basis(i,:)-basis(j,:))) == 1 && (A(i,j)==0)
            disp(strcat('Entry not right at ',num2str(i),',',num2str(j)));
        end
    end
end


end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A = generator5_circ(p0, q0, pr, qr, pl, ql, pc, qc)


% The order of the elements are as follows:
% 00000(1)
% 10000(2) 01000(3) 00100(4) 00010(5) 00001(6)
% 11000(7) 01100(8) 00110(9) 00011(10) 10001(11)
% 10100(12) 01010(13) 00101(14) 10010(15) 01001(16)
% 11100(17) 01110(18) 00111(19) 10011(20) 11001(21)
% 11010(22) 01101(23) 10110(24) 01011(25) 10101(26)
% 11110(27) 01111(28) 10111(29) 11011(30) 11101(31)
% 11111(32)



A(1:32,1:32)=0;
A(1,2:6)=p0;


A(2,1)=q0; A(2,7)=pr; A(2,11)=pl; A(2,12)=p0; A(2,15)=p0;
A(3,1)=q0; A(3,8)=pr; A(3,7)=pl; A(3,13)=p0; A(3,16)=p0;
A(4,1)=q0; A(4,9)=pr; A(4,8)=pl; A(4,14)=p0; A(4,12)=p0;
A(5,1)=q0; A(5,10)=pr; A(5,9)=pl; A(5,15)=p0; A(5,13)=p0;
A(6,1)=q0; A(6,11)=pr; A(6,10)=pl; A(6,16)=p0; A(6,14)=p0;

A(7,2)=qr; A(7,3)=ql; A(7,17)=pr; A(7,21)=pl; A(7,22)=p0;
A(8,3)=qr; A(8,4)=ql; A(8,18)=pr; A(8,17)=pl; A(8,23)=p0;
A(9,4)=qr; A(9,5)=ql; A(9,19)=pr; A(9,18)=pl; A(9,24)=p0;
A(10,5)=qr; A(10,6)=ql; A(10,20)=pr; A(10,19)=pl; A(10,25)=p0;
A(11,6)=qr; A(11,2)=ql; A(11,21)=pr; A(11,20)=pl; A(11,26)=p0;

A(12,2)=q0; A(12,4)=q0; A(12,17)=pc; A(12,24)=pr; A(12,26)=pl;
A(13,3)=q0; A(13,5)=q0; A(13,18)=pc; A(13,25)=pr; A(13,22)=pl;
A(14,4)=q0; A(14,6)=q0; A(14,19)=pc; A(14,26)=pr; A(14,23)=pl;
A(15,5)=q0; A(15,2)=q0; A(15,20)=pc; A(15,22)=pr; A(15,24)=pl;
A(16,6)=q0; A(16,3)=q0; A(16,21)=pc; A(16,23)=pr; A(16,25)=pl;


A(17,7)=qr; A(17,8)=ql; A(17,12)=qc; A(17,27)=pr; A(17,31)=pl;
A(18,8)=qr; A(18,9)=ql; A(18,13)=qc; A(18,28)=pr; A(18,27)=pl;
A(19,9)=qr; A(19,10)=ql; A(19,14)=qc; A(19,29)=pr; A(19,28)=pl;
A(20,10)=qr; A(20,11)=ql; A(20,15)=qc; A(20,30)=pr; A(20,29)=pl;
A(21,11)=qr; A(21,7)=ql; A(21,16)=qc; A(21,31)=pr; A(21,30)=pl;

A(22,7)=q0; A(22,15)=qr; A(22,13)=ql; A(22,27)=pc; A(22,30)=pc;
A(23,8)=q0; A(23,16)=qr; A(23,14)=ql; A(23,28)=pc; A(23,31)=pc;
A(24,9)=q0; A(24,12)=qr; A(24,15)=ql; A(24,29)=pc; A(24,27)=pc;
A(25,10)=q0; A(25,13)=qr; A(25,16)=ql; A(25,30)=pc; A(25,28)=pc;
A(26,11)=q0; A(26,14)=qr; A(26,12)=ql; A(26,31)=pc; A(26,29)=pc;

A(27,24)=qc; A(27,22)=qc; A(27,18)=ql; A(27,17)=qr; A(27,32)=pc;
A(28,25)=qc; A(28,23)=qc; A(28,19)=ql; A(28,18)=qr; A(28,32)=pc;
A(29,26)=qc; A(29,24)=qc; A(29,20)=ql; A(29,19)=qr; A(29,32)=pc;
A(30,22)=qc; A(30,25)=qc; A(30,21)=ql; A(30,20)=qr; A(30,32)=pc;
A(31,23)=qc; A(31,26)=qc; A(31,17)=ql; A(31,21)=qr; A(31,32)=pc;

A(32,27:31)=qc;

A = A - diag(sum(A,2)); % Make the rows sum to zero. Conservation of probability flow.



end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

