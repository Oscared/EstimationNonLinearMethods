%Generate process noise 
%Type 0: gaussian noise
%Type 1: small, random noise every 10th dt_sim
%Type 2: specific noise that is randomly distributed
function noise = gen_noise(th, type, aux)

sign_size = size(th);

if type==0
    if (nargin < 3)
        noise = sqrt(0.1).*randn(sign_size);
    else
        noise = sqrt(diag(aux))'.*randn(sign_size);
    end
elseif type==1
    noise = zeros(sign_size);
    for i=1:max(sign_size)
        if mod(i,10) == 0
            noise(i,:) = 0.5.*randn(size(noise(i,:)));
        end
    end
elseif type==2
    noise = zeros(sign_size);
    for i=1:max(sign_size)
        if rand > 0.5
            noise(i,:) = 0.7.*ones(size(noise(i,:)));
        end
    end
elseif type==3
    noise = zeros(sign_size);
    for i=1:max(sign_size)
        if th(i,1)>pi-0.1 || th(i,1)<pi+0.1
            th(i,2) = th(i,2) + 0.2;
        end
    end
end
end