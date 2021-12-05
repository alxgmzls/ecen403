

clear all;
close all;

syms K

K = 4;
K_Sq = K^2;
K_fact = factorial(K);

r = randi([1 10000],1,K);

sortedR = sort(r,'descend');

POP = sortedR.*(1/min(r))

% popularityVect = [2 1 1 1];

popularityVect = POP;

[PMF] = generatePMF_K(popularityVect,K);

% initialize marginal PMF for W
for i=1:height(PMF)
    
    marginalW(i) = sum(PMF(i,:));
    
end

[qType1, qType2] = genQueryStruct(K);

[fPattern1, fPattern2] = gen_fPatterns(K);
        
%----------------------------------------
allPermutations = perms([1:1:K]);

qTypes_fTypes = cell(2,2);

qTypes = {qType1; qType2};
fTypes = {fPattern1; fPattern2};

qTypes_fTypes = horzcat(qTypes,fTypes);

% there are 2*(K!)*(K^2) variables
numVariables = 2*factorial(K)*(K^2); % a 1 1 1 1

numVarsInQuery = numVariables/2;
%lowerbound of P values are all 0
lb = sparse(zeros(1,numVariables));
%upperbound of P values initially set to 1, but some will be set to 0
ub = ones(1,numVariables);

% this serves as a reminder of the convention
% w and s according to relative positions of an index 1 to K^2
% relative to subsets of K_Sq that characterize a distinct Query
% Example for K = 5
    % w,s = 1,1     has relative position 1
    % w,s = 1,2     has relative position 2
    % w,s = 1,3     has relative position 3
    % w,s = 1,4     has relative position 4
    % w,s = 1,5     has relative position 5
    
    % w,s = 2,1     has relative position 6
    % w,s = 2,2     has relative position 7
    % w,s = 2,3     has relative position 8
    % w,s = 2,4     has relative position 9
    % w,s = 2,5     has relative position 10
    
    % w,s = 3,1     has relative position 11
    % w,s = 3,2     has relative position 12
    % w,s = 3,3     has relative position 13
    % w,s = 3,4     has relative position 14
    % w,s = 3,5     has relative position 15
    
    % w,s = 4,1     has relative position 16
    % w,s = 4,2     has relative position 17
    % w,s = 4,3     has relative position 18
    % w,s = 4,4     has relative position 19
    % w,s = 4,5     has relative position 20
    
    % w,s = 5,1     has relative position 21
    % w,s = 5,2     has relative position 22
    % w,s = 5,3     has relative position 23
    % w,s = 5,4     has relative position 24
    % w,s = 5,5     has relative position 25

%loop through each Px, and set to 0 infeasible combinations of W and S
for x=1:width(ub)
    
    %first, extract query type
    QTYPE_INDEX = ceil(x/((K_Sq)*K_fact));
    
    %then, extract feasability pattern
    FPATTERN = qTypes_fTypes{QTYPE_INDEX,2};
    
    %now, we must extract w and s from our convention
    %first, "normalize" index into size 25
    rel_pos = mod(x-1,K_Sq)+1;
    
    %from, relative position, find out w and s candidates
    W = ceil(rel_pos/K);
    S = mod(rel_pos-1,K)+1;
    
    
    %we need the permutation index
    PERM_INDEX = mod(ceil(x/K_Sq)-1,K_fact)+1;
    PERM = allPermutations(PERM_INDEX,:);
    
    
    if(QTYPE_INDEX == 2 && PERM_INDEX ~= K_fact)
       ub(1,x) = 0;
       continue;
    end
    
%     if(ismember(PERM_INDEX,UsePermutations)==false)
%         ub(1,x) = 0;
%         continue;
%     end
%     
    %now we have everything to call isFeasible
    if(isFeasible(FPATTERN,PERM,W,S) == 0)
        ub(1,x) = 0; %if non-feasible, set element in indexed column to 0
    end
%     else
%         fprintf('FPATTERN %i, PERM %s, W,S = %i,%i isFeasible = %i \n', QTYPE_INDEX,num2str(PERM),W,S,isFeasible(FPATTERN,PERM,W,S));
%     end
    
end

% defining Aeq
[base] = populateBaseBlock(zeros(K,K_Sq), marginalW, PMF, K_Sq, K);
sparse_base = sparse(base);
clear base;

numPijs = numVariables/K_Sq;
Ac = repmat({sparse_base}, 1, numPijs);  
Aeq = sparse(blkdiag(Ac{:}));

%-----------------------------------------------------------
%-----------------------------------------------------------

AeqAppend = [];

for j=1:K_Sq
    temp = zeros(1,numVariables);
    for i=1:width(temp)
        if(mod(i-1,K_Sq)+1 == j)
            temp(1,i) = 1;
        end
    end
    AeqAppend = [AeqAppend;temp];
end

AeqAppend = sparse(AeqAppend);

%----------------------------------------------
%------------ putting them together -------------
%----------------------------------------------

Aeq = vertcat(Aeq,AeqAppend);

%now we must develop optimization condition, or we want to minimize this
%expression. We will use the last KxK rows of Aeq, or the rows due to
%condition Two.

downloadOptimization = zeros(1,numVariables);
conditionTwo = AeqAppend;
clear AeqAppend

for row=1:height(conditionTwo)
    % for each row, scale GPC code by number of rows in qType1, and number
    % of rows in qType2
    for col=1:width(conditionTwo)
        if(col <= numVarsInQuery)
            conditionTwo(row,col) = conditionTwo(row,col)*height(qType1);
        else
            conditionTwo(row,col) = conditionTwo(row,col)*height(qType2);
        end
    end
end

conditionTwo = sparse(conditionTwo);

%num2str(conditionTwo)

for row=1:height(conditionTwo)
    %extract relative position from convention and column index
    rel_pos = mod(row-1,K_Sq)+1;
    pmfW_index = ceil(rel_pos/K);
    pmfS_index = mod(rel_pos-1,K)+1;
        
    %scale row by PMF
    conditionTwo(row,:) = conditionTwo(row,:)*PMF(pmfW_index,pmfS_index);
end

%num2str(conditionTwo)

for index=1:width(downloadOptimization)
    downloadOptimization(1,index) = sum(conditionTwo(:,index));
end


% need to create matrix of form 0 1 1 1 ... 1
%                              1 0 1 1 ... 1
%                              1 1 0 1 ... 1
%                              1 1 1 0 ... 1
% or an all ones matrix with zero diagonal

K_K_ones = ones(K,K);
normaliz_constr = K_K_ones - diag(diag(K_K_ones));

numPrivacyCondEqns = K*numPijs;
Beq = sparse(vertcat(zeros(numPrivacyCondEqns,1),normaliz_constr(:)));

%now we call linprog, storing the values of Px in X, and the # of downloads
%in Z


[X,Z,exitFlag] = linprog(downloadOptimization,[],[],Aeq,Beq,lb,ub);

NZp = [];
for index=1:height(X)
    if(X(index,1) ~= 0)
        NZp = [NZp;X(index,1) index mod(index-1,K_Sq)+1];
    end
end

NZp =sortrows(NZp,3);

PsUsed = [];

for index=1:height(NZp)
    rel_pos = mod(NZp(index,2)-1,K_Sq)+1;
    W = ceil(rel_pos/K);
    S = mod(rel_pos-1,K)+1;
    Q = ceil(NZp(index,2)/(numVarsInQuery));
    P = mod(ceil(NZp(index,2)/K_Sq)-1,K_fact)+1;
    
    if(ismember(P,PsUsed)==false)
       PsUsed = [PsUsed,P]; 
    end
    
%     if(Q == 1)
%         probability2row = probability2row + PWS(W,S)*NZp(index,1);
%     else
%         probability3row = probability3row + PWS(W,S)*NZp(index,1);
%     end
    
%     if(Q == 1)
%         exact_sendQuery1WithProb_matrix(W,S) = NZp(index,1);
%         normalized_pij(W,S) = (NZp(index,1)/((marginalW(1,W)/marginalW(1,4))*(0.25/PWS(W,S))))*(POP(1)+POP(2)+POP(3));
%         
%         if((W == 1 && S == 2) | (W == 2 && S == 1) | (W == 4 && S == 3) | (W == 3 && S == 4))
%             
%                 normalized_pij(W,S) = normalized_pij(W,S)*((marginalW(1,3)*POP(2))/(POP(3)*marginalW(1,2)))*(marginalW(1,2)/(POP(2)*marginalW(1,4)));
%                 
%         elseif((W == 1 && S == 3) | (W == 3 && S == 1) | (W == 2 && S == 4) | (W == 4 && S == 2))
%                 
%                 %normalized_pij(W,S) = normalized_pij(W,S)*((marginalW(1,3)*POP(4))/(POP(3)*marginalW(1,4)))*(marginalW(1,2)/(POP(2)*marginalW(1,4)));
%                 %normalized_pij(W,S) = normalized_pij(W,S)*(marginalW(1,3)/(POP(3)*marginalW(1,4)))*((marginalW(1,2)*POP(3))/(POP(2)*marginalW(1,3)));
%                 normalized_pij(W,S) = normalized_pij(W,S)*(marginalW(1,2)/(POP(2)*marginalW(1,4)));
%             
%         elseif((W == 4 && S == 1) | (W == 3 && S == 2) | (W == 2 && S == 3) | (W == 1 && S == 4))
%             % both of these currently suck
%                   option1 = normalized_pij(W,S)*(marginalW(1,2)/(POP(2)*marginalW(1,4)))*((POP(1)+POP(2)+POP(4))/(POP(1)+POP(2)+POP(3)));
%                   %option2 = normalized_pij(W,S)*(marginalW(1,2)/(POP(2)*marginalW(1,4)));
%                   option3 = normalized_pij(W,S)*(marginalW(1,2)/(POP(2)*marginalW(1,4)))*((POP(2)*marginalW(1,1))/(POP(1)*marginalW(1,2)));
%                   
%                   normalized_pij(W,S) = max(option1,option3);
% 
% 
%         end
%     end
    %fprintf('for W,S = %d,%d , the probability of choosing Q-type %d should be %s. \n',W,S,Q,rats(NZp(index,1)));
    fprintf('P(QTYPE %i, PERM %i, W = %i, S = %i) = %s \n',Q,P,W,S,rats(NZp(index,1)));
    %disp(allPermutations(P,:))
    %allPermutations(P,:)
end

fprintf('Expected download cost: %f \n',Z);
% rats(Z)

rats(PMF)
% 
% for i=1:width(PsUsed)
%     
%     disp(allPermutations(PsUsed(i),:));
%     
% end


function [PWS] = generatePMF_K(POP,K)

    PWS = zeros(K,K);
    
    for i=1:height(PWS)
        
        for j=1:width(PWS)
            
            if(i~=j)
                
                PWS(i,j) = (1/K)*(POP(i)/(sum(POP)-POP(j)));
                
            end
            
        end
        
    end

end

function [Aeq_base] = populateBaseBlock(empty, marginalW, PWS, K_Sq, K)

    Aeq_base = empty;

    for i=1:height(Aeq_base)
       
        for j=1:width(Aeq_base)
            
            rel_pos = mod(j-1,K_Sq)+1;
            
            PWS_i = ceil(rel_pos/K);
            PWS_j = mod(rel_pos-1,K)+1;
            
            if(ceil(rel_pos/K) == i)    
                
                Aeq_base(i,j) = PWS(PWS_i,PWS_j)*(1-marginalW(i));
                
            else
                
                Aeq_base(i,j) = -1*PWS(PWS_i,PWS_j)*marginalW(i);
                
            end
            
        end
        
    end
    
end

function [feasability] = isFeasible(feasabilityPattern,qPerm,w,s)
    

    temp = feasabilityPattern;
    
    %scale columns by element in permutation
    for col=1: width(qPerm);
        temp(:,col) = feasabilityPattern(:,col)*qPerm(1,col);
    end
    
    %add transposed perm column horizontally
    
    temp = horzcat(transpose(qPerm),temp);
    
    %identify which row (from first column of temp) contains the candidate
    %demand
    
    si_row = find(temp(:,1) == w);
    
    %check if elements in row (without the first column) contain candidate
    %SI
    
    if(ismember(s,temp(si_row,2:end))==true)
        feasability = 1; % yes, the candidate demand and candidate side information pairs are valid
    else
        feasability = 0; % no, the pairing is non-feasible
    end

    
end

function [GPC, MDS] = genQueryStruct(K)

    GPC = [];

    for d = 1:floor(K/2);

        GPC = blkdiag(GPC, [1 1]);

    end

    if(mod(K,2)==1)

        GPC = blkdiag(GPC, [1]);

    end

    MDS = horzcat(eye(K-1),ones(K-1,1));

end

function [fP1, fP2] = gen_fPatterns(K)

    fP1 = [];
    
    fP2 = [];
    
    baseBlock = [0 1;
                 1 0];
    
    for i=1:floor(K/2)
        
        fP1 = blkdiag(fP1,baseBlock);
        
    end
   
    if(mod(K,2) == 1)
        
        fP1 = [fP1 zeros(K-1,1)];
        lastVect = ones(1,K);
        lastVect(1,K) = 0;
        fP1 = vertcat(fP1,lastVect);
        
    end
    
    fP2 = ones(K,K);
    fP2 = fP2 - diag(diag(fP2));
    
    
end

