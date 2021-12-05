close all;
clear all;

syms a k

data = [];

for k = 4:300

    a = 10000;

    popularityVect = [a ones(1,k-1)];

    if mod(k,2) == 0 % if k is even

        syms q p r

        x = k - 1;

        y = a + k - 2;

        param1 = x/y;

        param2 = ((1/x) + (k-2)/y);

        r = q * (1/param1);

        p = q * (1/param2);

        cost = (k-1)*(a/y)*(1/k)*(p*(k/2)+(1-p)*(k-1))...
             + (k-1)*(1/x)*(1/k)*(q*(k/2)+(1-q)*(k-1))...
     + 2*nchoosek(k-1,2)*(1/y)*(1/k)*(r*(k/2)+(1-r)*(k-1));

        cost_coeffs = coeffs(cost);

        bias = double(cost_coeffs(1));

        savings_term = double(cost_coeffs(2));

        lb = [0];

        ub = [min([param1 param2 1])];

        [q,fval] = linprog([savings_term],[],[],[],[],lb,ub);

        r = q * (1/param1);

        p = q * (1/param2);

         expected_cost = fval + bias;
%         if (k <= 10)
%             disp(popularityVect);
%         end
%         fprintf('for W = (a), S = (1), p = %s \n',rats(p));
%         fprintf('for W = (1), S = (a), q = %s \n',rats(q));
%         fprintf('for W = (1), S = (1), r = %s \n',rats(r));
%         fprintf('Expected Download Cost: %f \n', expected_cost);
%         fprintf('Cost Savings: %f \n',fval*-1);

          percentCostSavings = ((k-1)-expected_cost)/(k-1);
          
          

    else % else k is odd

        syms q1 q2 p1 p2 r1 r2

        x = k - 1;

        y = a + k - 2;

        param1 = ((1/x) + (k-2)/y);

        param2 = (param1*x)/(k-2);

        param3 = (x/y);

        param4 = (param1)/(1 + (k-2)/(k-1));

        p1 = q1 * (1/param1);

        p2 = q1 * (1/param2);

        r1 = q1 * (1/param3);

        r2 = r1 - q2 * (y/(x*(k-2)));

        param5 = (x/(2*y)) + q2/(2*(k-2));

    %     lb = [max(q2/(k-2),0) 0];
        lb = [0 0];

    %     ub = [min([param1, param2, param3, param4, param5, q2, 1])  1];

        ub = [min([param1, param2, param3, param4, 1])  1];

        cost = (k-1)*(a/y)*(1/k)*((p1+p2)*((floor(k/2))+1)+(1-p1-p2)*(k-1))...
                        +  (1/k)*((q1+q2)*((floor(k/2))+1)+(1-q1-q2)*(k-1))...
    + 2*nchoosek(k-1,2)*(1/y)*(1/k)*((r1+r2)*((floor(k/2))+1)+(1-r1-r2)*(k-1));

        cost_coeffs = coeffs(cost);

        bias = double(cost_coeffs(1));

        if(width(cost_coeffs) <= 2 )

            savings_term = double(cost_coeffs(2));

            [B,fval] = linprog([savings_term 0],[1 1; -1 -1; (y/x) -(y/(x*(k-2))); 1 -1; 0 1;-(y/x) (y/(x*(k-2))); -(k-2) 1; k-2 -1],[1 0 1 0 1 0 0 (x*(k-2))/(y)],[],[],lb,ub);

        else

            savings_term2 = double(cost_coeffs(2));

            savings_term1 = double(cost_coeffs(3));

            [B,fval] = linprog([savings_term1 savings_term2],[1 1; -1 -1; (y/x) -(y/(x*(k-2))); 1 -1; 0 1;-(y/x) (y/(x*(k-2))); -(k-2) 1; k-2 -1],[1 0 1 0 1 0 0 (x*(k-2))/(y)],[],[],lb,ub);


        end


    %     [B,fval] = linprog([savings_term 0],[-(k-2) 1; k-2 -1; 1 1 ; -1 -1],[0 (x*(k-2))/(y) 1 0],[],[],lb,ub);

        q1 = B(1);

        q2 = B(2);

        expected_cost = fval + bias;

        p1 = q1 * (1/param1);

        p2 = q1 * (1/param2);

        r1 = q1 * (1/param3);

        r2 = r1 - q2 * (y/(x*(k-2)));

        cost = (k-1)*(a/y)*(1/k)*((p1+p2)*((floor(k/2))+1)+(1-p1-p2)*(k-1))...
                        +  (1/k)*((q1+q2)*((floor(k/2))+1)+(1-q1-q2)*(k-1))...
    + 2*nchoosek(k-1,2)*(1/y)*(1/k)*((r1+r2)*((floor(k/2))+1)+(1-r1-r2)*(k-1));


%         if (k <= 10)
%             disp(popularityVect);
%         end
%         fprintf('for W = (a), S = (1), p1 = %s \n',rats(p1));
%         fprintf('for W = (a), S = (1), p2 = %s \n',rats(p2));
%         fprintf('for W = (1), S = (a), q1 = %s \n',rats(q1));
%         fprintf('for W = (1), S = (a), q2 = %s \n',rats(q2));
%         fprintf('for W = (1), S = (1), r1 = %s \n',rats(r1));
%         fprintf('for W = (1), S = (1), r2 = %s \n',rats(r2));
%         fprintf('Expected Download Cost: %f \n', cost);
%         fprintf('Cost Savings: %f \n',(k-1)-cost);
    
          percentCostSavings = ((k-1)-cost)/(k-1);
    end
    
    
    data = [data; k percentCostSavings];
    
end

scatter(data(:,1),data(:,2))