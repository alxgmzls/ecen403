

for test_i = 0:100
    
    
    
    r = randi([1 100],1,4);

    sortedR = sort(r,'descend');

    POP = sortedR.*(1/min(r))



    % POP = [102 1 1 1]

    UseQueries = [1,5];


    UsePermutations = [20 22 24];

    PWS = zeros(4,4);
    PWS(1,:) = [0 POP(1)/sum(POP([1,3,4])) POP(1)/sum(POP([1,2,4])) POP(1)/sum(POP([1,2,3]))]/4;
    PWS(2,:) = [POP(2)/sum(POP([2,3,4])) 0 POP(2)/sum(POP([1,2,4])) POP(2)/sum(POP([1,2,3]))]/4;
    PWS(3,:) = [POP(3)/sum(POP([2,3,4])) POP(3)/sum(POP([1,3,4])) 0 POP(3)/sum(POP([1,2,3]))]/4;
    PWS(4,:) = [POP(4)/sum(POP([2,3,4])) POP(4)/sum(POP([1,3,4])) POP(4)/sum(POP([1,2,4])) 0]/4;

    POPS = zeros(4,4);
    POPS(1,:) = [0 POP(1)/sum(POP([1,3,4])) POP(1)/sum(POP([1,2,4])) POP(1)/sum(POP([1,2,3]))];
    POPS(2,:) = [POP(2)/sum(POP([2,3,4])) 0 POP(2)/sum(POP([1,2,4])) POP(2)/sum(POP([1,2,3]))];
    POPS(3,:) = [POP(3)/sum(POP([2,3,4])) POP(3)/sum(POP([1,3,4])) 0 POP(3)/sum(POP([1,2,3]))];
    POPS(4,:) = [POP(4)/sum(POP([2,3,4])) POP(4)/sum(POP([1,3,4])) POP(4)/sum(POP([1,2,4])) 0];
    % all "types" of query structures

    qType1 = [1 1 0 0;
              0 0 1 1];

    qType2 = [1 0 0 0;
              0 1 1 0; 
              0 0 0 1];

    qType3 = [1 0 0 0; 
              0 1 1 0;
              0 0 1 1];

    qType4 = [1 0 0 0;
              1 1 0 0;
              0 0 1 1];

    qType5 = [1 1 0 0;
              0 1 1 0;
              0 0 1 1];

    qType6 = [1 1 0 0;
              1 0 1 0;
              1 0 0 1];

    qType7 = [1 0 0 0;
              0 1 0 0;
              0 0 1 0;
              0 0 0 1];

    % all recoverability patterns for same qType

    fPattern1 = [0 1 0 0;
                 1 0 0 0;
                 0 0 0 1;
                 0 0 1 0];

    fPattern2 = [0 1 1 1;
                 0 0 1 0;
                 0 1 0 0;
                 1 1 1 0];

    fPattern3 = [0 1 1 1;
                 0 0 1 1;
                 0 1 0 1;
                 0 1 1 0];

    fPattern4 = [0 1 1 1;
                 1 0 1 1;
                 0 0 0 1;
                 0 0 1 0];

    fPattern5 = [0 1 1 1;
                 1 0 1 1;
                 1 1 0 1;
                 1 1 1 0];

    fPattern6 = [0 1 1 1;
                 1 0 1 1;
                 1 1 0 1;
                 1 1 1 0];

    fPattern7 = [0 1 1 1;
                 1 0 1 1;
                 1 1 0 1;
                 1 1 1 0];

    % all permutations of possible message placements
    v = [1 2 3 4];
    allPermutations = perms(v);

    qTypes_fTypes = cell(7,2);

    qTypes = {qType1; qType2; qType3; qType4; qType5; qType6; qType7};
    fTypes = {fPattern1; fPattern2; fPattern3; fPattern4; fPattern5; fPattern6; fPattern7};

    qTypes_fTypes = horzcat(qTypes,fTypes);

    %lowerbound of P values are all 0
    lb = zeros(1,2688);
    %upperbound of P values initially set to 1, but some will be set to 0
    ub = ones(1,2688);

    % this serves as a reminder of the convention
    % w and s according to relative positions of an index 1 to 2688
    % relative to subsets of 16 that characterize a distinct Query

        % w,s = 1,1     has relative position 1
        % w,s = 1,2     has relative position 2
        % w,s = 1,3     has relative position 3
        % w,s = 1,4     has relative position 4

        % w,s = 2,1     has relative position 5
        % w,s = 2,2     has relative position 6
        % w,s = 2,3     has relative position 7
        % w,s = 2,4     has relative position 8

        % w,s = 3,1     has relative position 9
        % w,s = 3,2     has relative position 10
        % w,s = 3,3     has relative position 11
        % w,s = 3,4     has relative position 12

        % w,s = 4,1     has relative position 13
        % w,s = 4,2     has relative position 14
        % w,s = 4,3     has relative position 15
        % w,s = 4,4     has relative position 16



    %loop through each Px, and set to 0 infeasible combinations of W and S
    for x=1:width(ub)

        %first, extract query type
        QTYPE_INDEX = ceil(x/(16*24));
        % QTYPE = qTypes_fTypes{QTYPE_INDEX,1}; we dont really need QTYPE

        %if we dont want to use the QTYPE, we can restrict ourselves by upper
        %bounding the corresponding variables to 0 and continuing
        if(ismember(QTYPE_INDEX,UseQueries)==false)
            ub(1,x) = 0;
            continue;
        end

        %then, extract feasability pattern
        FPATTERN = qTypes_fTypes{QTYPE_INDEX,2};

        %now, we must extract w and s from our convention
        %first, "normalize" index into size 16 (1,2..15,16,1,2...15,16...)
        rel_pos = mod(x-1,16)+1;

        %from, relative position, find out w and s candidates
        W = ceil(rel_pos/4);
        S = mod(rel_pos-1,4)+1;

    %     if(ismember(S,UseSI) == false | ismember(W,UseWant) == false)
    %         ub(1,x) = 0;
    %         continue;
    %     end

        %we need the permutation index
        % ceil(x/16) divides Px's into 168 subsets, which each correspond to
        % distinct query. Amongst those 168 subsets, taking sequential subsets
        % of size 24 "lumps" same query types, since there are 1,2...24
        % permutations for K=4,
        PERM_INDEX = mod(ceil(x/16)-1,24)+1;
        PERM = allPermutations(PERM_INDEX,:);

        if(QTYPE_INDEX == 5 && PERM_INDEX ~= 24)
           ub(1,x) = 0;
           continue;
        end

        if(ismember(PERM_INDEX,UsePermutations)==false)
            ub(1,x) = 0;
            continue;
        end

        %now we have everything to call isFeasible
        if(isFeasible(FPATTERN,PERM,W,S) == 0)
            ub(1,x) = 0; %if non-feasible, set element in indexed column to 0
        end
    %     else
    %         fprintf('FPATTERN %i, PERM %s, W,S = %i,%i isFeasible = %i \n', QTYPE_INDEX,num2str(PERM),W,S,isFeasible(FPATTERN,PERM,W,S))
    %     end

    end
    % defining Aeq

    base = zeros(4,16);


    % w = 1 equations
    base(1,1)  = PWS(1,1)*(1-sum(PWS(1,:)));
    base(1,2)  = PWS(1,2)*(1-sum(PWS(1,:)));
    base(1,3)  = PWS(1,3)*(1-sum(PWS(1,:)));
    base(1,4)  = PWS(1,4)*(1-sum(PWS(1,:)));

    base(1,5)  = PWS(2,1)*sum(PWS(1,:))*-1;
    base(1,6)  = PWS(2,2)*sum(PWS(1,:))*-1;
    base(1,7)  = PWS(2,3)*sum(PWS(1,:))*-1;
    base(1,8)  = PWS(2,4)*sum(PWS(1,:))*-1;

    base(1,9)  = PWS(3,1)*sum(PWS(1,:))*-1;
    base(1,10) = PWS(3,2)*sum(PWS(1,:))*-1;
    base(1,11) = PWS(3,3)*sum(PWS(1,:))*-1;
    base(1,12) = PWS(3,4)*sum(PWS(1,:))*-1;

    base(1,13) = PWS(4,1)*sum(PWS(1,:))*-1;
    base(1,14) = PWS(4,2)*sum(PWS(1,:))*-1;
    base(1,15) = PWS(4,3)*sum(PWS(1,:))*-1;
    base(1,16) = PWS(4,4)*sum(PWS(1,:))*-1;

    %-------------------------------------------------------
    %-------------------------------------------------------

    % w = 2 equations

    base(2,1)  = PWS(1,1)*sum(PWS(2,:))*-1;
    base(2,2)  = PWS(1,2)*sum(PWS(2,:))*-1;
    base(2,3)  = PWS(1,3)*sum(PWS(2,:))*-1;
    base(2,4)  = PWS(1,4)*sum(PWS(2,:))*-1;

    base(2,5)  = PWS(2,1)*(1-sum(PWS(2,:)));
    base(2,6)  = PWS(2,2)*(1-sum(PWS(2,:)));
    base(2,7)  = PWS(2,3)*(1-sum(PWS(2,:)));
    base(2,8)  = PWS(2,4)*(1-sum(PWS(2,:)));

    base(2,9)  = PWS(3,1)*sum(PWS(2,:))*-1;
    base(2,10) = PWS(3,2)*sum(PWS(2,:))*-1;
    base(2,11) = PWS(3,3)*sum(PWS(2,:))*-1;
    base(2,12) = PWS(3,4)*sum(PWS(2,:))*-1;

    base(2,13) = PWS(4,1)*sum(PWS(2,:))*-1;
    base(2,14) = PWS(4,2)*sum(PWS(2,:))*-1;
    base(2,15) = PWS(4,3)*sum(PWS(2,:))*-1;
    base(2,16) = PWS(4,4)*sum(PWS(2,:))*-1;

    %-------------------------------------------------------
    %-------------------------------------------------------

    % w = 3 equations

    base(3,1)  = PWS(1,1)*sum(PWS(3,:))*-1;
    base(3,2)  = PWS(1,2)*sum(PWS(3,:))*-1;
    base(3,3)  = PWS(1,3)*sum(PWS(3,:))*-1;
    base(3,4)  = PWS(1,4)*sum(PWS(3,:))*-1;

    base(3,5)  = PWS(2,1)*sum(PWS(3,:))*-1;
    base(3,6)  = PWS(2,2)*sum(PWS(3,:))*-1;
    base(3,7)  = PWS(2,3)*sum(PWS(3,:))*-1;
    base(3,8)  = PWS(2,4)*sum(PWS(3,:))*-1;

    base(3,9)  = PWS(3,1)*(1-sum(PWS(3,:)));
    base(3,10) = PWS(3,2)*(1-sum(PWS(3,:)));
    base(3,11) = PWS(3,3)*(1-sum(PWS(3,:)));
    base(3,12) = PWS(3,4)*(1-sum(PWS(3,:)));

    base(3,13) = PWS(4,1)*sum(PWS(3,:))*-1;
    base(3,14) = PWS(4,2)*sum(PWS(3,:))*-1;
    base(3,15) = PWS(4,3)*sum(PWS(3,:))*-1;
    base(3,16) = PWS(4,4)*sum(PWS(3,:))*-1;

    %-------------------------------------------------------
    %-------------------------------------------------------

    % w = 4

    base(4,1)  = PWS(1,1)*sum(PWS(4,:))*-1;
    base(4,2)  = PWS(1,2)*sum(PWS(4,:))*-1;
    base(4,3)  = PWS(1,3)*sum(PWS(4,:))*-1;
    base(4,4)  = PWS(1,4)*sum(PWS(4,:))*-1;

    base(4,5)  = PWS(2,1)*sum(PWS(4,:))*-1;
    base(4,6)  = PWS(2,2)*sum(PWS(4,:))*-1;
    base(4,7)  = PWS(2,3)*sum(PWS(4,:))*-1;
    base(4,8)  = PWS(2,4)*sum(PWS(4,:))*-1;

    base(4,9)  = PWS(3,1)*sum(PWS(4,:))*-1;
    base(4,10) = PWS(3,2)*sum(PWS(4,:))*-1;
    base(4,11) = PWS(3,3)*sum(PWS(4,:))*-1;
    base(4,12) = PWS(3,4)*sum(PWS(4,:))*-1;

    base(4,13) = PWS(4,1)*(1-sum(PWS(4,:)));
    base(4,14) = PWS(4,2)*(1-sum(PWS(4,:)));
    base(4,15) = PWS(4,3)*(1-sum(PWS(4,:)));
    base(4,16) = PWS(4,4)*(1-sum(PWS(4,:)));

    %-------------------------------------------------------
    %-------------------------------------------------------

    Ac = repmat({base}, 1, 168);  
    Aeq = blkdiag(Ac{:});

    %-----------------------------------------------------------
    %-----------------------------------------------------------

    AeqAppend = [];

    for j=1:16
        temp = zeros(1,2688);
        for i=1:width(temp)
            if(mod(i-1,16)+1 == j)
                temp(1,i) = 1;
            end
        end
        AeqAppend = [AeqAppend;temp];
    end


    %----------------------------------------------
    %------------ putting them together -------------
    %----------------------------------------------

    Aeq = vertcat(Aeq,AeqAppend);

    %now we must develop optimization condition, or we want to minimize this
    %expression. We will use the last KxK rows of Aeq, or the rows due to
    %condition Two.
    downloadOptimization = zeros(1,2688);
    conditionTwo = AeqAppend;

    for row=1:height(conditionTwo)
        % for each row, scale columns related to 2, 3, and 4 downloads
        for col=1:width(conditionTwo)
            if(col <= 384)
                conditionTwo(row,col) = conditionTwo(row,col)*2;
            elseif(col > 384 && col <= 2304)
                conditionTwo(row,col) = conditionTwo(row,col)*3;
            elseif(col > 2304 && col <= 2688)
                conditionTwo(row,col) = conditionTwo(row,col)*4;
            end
        end
    end

    %num2str(conditionTwo)

    for row=1:height(conditionTwo)
        %extract relative position from convention and column index
        rel_pos = mod(row-1,16)+1;
        pmfW_index = ceil(rel_pos/4);
        pmfS_index = mod(rel_pos-1,4)+1;

        %scale row by PMF
        conditionTwo(row,:) = conditionTwo(row,:)*PWS(pmfW_index,pmfS_index);
    end

    %num2str(conditionTwo)

    for index=1:width(downloadOptimization)
        downloadOptimization(1,index) = sum(conditionTwo(:,index));
    end

    %num2str(downloadOptimization)

    %finally, we must construct the Beq matrix, which consists of a column
    %vector of 672 0's and 16 ones
    Beq = vertcat(zeros(672,1),[0;1;1;1
                                1;0;1;1;
                                1;1;0;1;
                                1;1;1;0]);

    %now we call linprog, storing the values of Px in X, and the # of downloads
    %in Z


    [X,Z,exitFlag] = linprog(downloadOptimization,[],[],Aeq,Beq,lb,ub);

    %NZp gathers non-zero P variables and row-lumps them with their index
    % NZp = [];
    % for index=1:height(X)
    %     if(X(index,1) ~= 0)
    %         NZp = [NZp;X(index,1) index];
    %     end
    % end

    marginalW = [];

    for i = 1:height(PWS)
        marginalW = [marginalW sum(PWS(i,:))];
    end

    NZp = [];
    for index=1:height(X)
        if(X(index,1) ~= 0)
            NZp = [NZp;X(index,1) index mod(index-1,16)+1];
        end
    end

    % this module orders NZp in such a way that prints "nicely"
    % added 7/30

    NZp =sortrows(NZp,3);

    probability2row = 0;
    probability3row = 0;

    PsUsed = [];

    exact_sendQuery1WithProb_matrix = zeros(4,4);
    normalized_pij = zeros(4,4);
    record = [];

    for index=1:height(NZp)
        rel_pos = mod(NZp(index,2)-1,16)+1;
        W = ceil(rel_pos/4);
        S = mod(rel_pos-1,4)+1;
        Q = ceil(NZp(index,2)/384);
        P = mod(ceil(NZp(index,2)/16)-1,24)+1;

        if(ismember(P,PsUsed)==false)
           PsUsed = [PsUsed,P]; 
        end

        if(Q == 1)
            probability2row = probability2row + PWS(W,S)*NZp(index,1);
        else
            probability3row = probability3row + PWS(W,S)*NZp(index,1);
        end

        if(Q == 1)
            exact_sendQuery1WithProb_matrix(W,S) = NZp(index,1);
            normalized_pij(W,S) = (NZp(index,1)/((marginalW(1,W)/marginalW(1,4))*(0.25/PWS(W,S))))*(POP(1)+POP(2)+POP(3));

            if((W == 1 && S == 2) | (W == 2 && S == 1) | (W == 4 && S == 3) | (W == 3 && S == 4))

                    normalized_pij(W,S) = normalized_pij(W,S)*((marginalW(1,3)*POP(2))/(POP(3)*marginalW(1,2)))*(marginalW(1,2)/(POP(2)*marginalW(1,4)));

            elseif((W == 1 && S == 3) | (W == 3 && S == 1) | (W == 2 && S == 4) | (W == 4 && S == 2))

                    %normalized_pij(W,S) = normalized_pij(W,S)*((marginalW(1,3)*POP(4))/(POP(3)*marginalW(1,4)))*(marginalW(1,2)/(POP(2)*marginalW(1,4)));
                    %normalized_pij(W,S) = normalized_pij(W,S)*(marginalW(1,3)/(POP(3)*marginalW(1,4)))*((marginalW(1,2)*POP(3))/(POP(2)*marginalW(1,3)));
                    normalized_pij(W,S) = normalized_pij(W,S)*(marginalW(1,2)/(POP(2)*marginalW(1,4)));

            elseif((W == 4 && S == 1) | (W == 3 && S == 2) | (W == 2 && S == 3) | (W == 1 && S == 4))
                % both of these currently suck
                      option1 = normalized_pij(W,S)*(marginalW(1,2)/(POP(2)*marginalW(1,4)))*((POP(1)+POP(2)+POP(4))/(POP(1)+POP(2)+POP(3)));
                      %option2 = normalized_pij(W,S)*(marginalW(1,2)/(POP(2)*marginalW(1,4)));
                      option3 = normalized_pij(W,S)*(marginalW(1,2)/(POP(2)*marginalW(1,4)))*((POP(2)*marginalW(1,1))/(POP(1)*marginalW(1,2)));

                      normalized_pij(W,S) = max(option1,option3);


            end
        end
        %fprintf('for W,S = %d,%d , the probability of choosing Q-type %d should be %s. \n',W,S,Q,rats(NZp(index,1)));
%        fprintf('P(QTYPE %i, PERM %i, W = %i, S = %i) = %s \n',Q,P,W,S,rats(NZp(index,1)));
        %allPermutations(P,:)
    end


%     fprintf('Expected download cost: %f \n',Z);
%     rats(Z)
% 
%     rats(PWS)

    approx_sendQuery1WithProb_matrix = zeros(4,4);


for i = 1:height(approx_sendQuery1WithProb_matrix)
    for j = 1:width(approx_sendQuery1WithProb_matrix)
        
        if( i~=j)
            
            indices = [1 2 3 4];
            
            k = min(setdiff(indices,[i j]));
            l = max(setdiff(indices,[i j]));
            
            if(i > j)
                
                r1 = PWS(j,i)/marginalW(j);
                r2 = PWS(k,l)/marginalW(k);
                
            elseif(i < j)
                    
                r1 = PWS(i,j)/marginalW(i);
                r2 = PWS(k,l)/marginalW(k);
                
            end
            
            approx_sendQuery1WithProb_matrix(i,j) = min(r1, r2)*(marginalW(i)/PWS(i,j));
            
            pij = min(r1, r2)*(marginalW(i)/PWS(i,j));
            
            if(pij > 1.00001 || pij < 0)
                fprintf('invalid pij...stopping program... \n');
                return;
            end
        end
    end
end

    % computation = (1/marginalW(1,4))*(marginalW(1,1)*0.25)*(1/PWS(1,2))*(1/(POP(1,1)+POP(1,2)+POP(1,3)));
    % rats(computation)

%     for i = 1:width(marginalW)
%         fprintf('probability that W = %i : %s\n',i, rats(marginalW(1,i)));
%     end

    if(abs(approx_sendQuery1WithProb_matrix-exact_sendQuery1WithProb_matrix)<0.0001)
        %fprintf('approximation is a perfect fit for this lambda \n');
    else
        fprintf('approximation is flawed for this lambda \n');  
        fprintf('...stopping program... \n');
        disp(POP);
    end

    %rats(NZp)
    %rats(normalized_pij)
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