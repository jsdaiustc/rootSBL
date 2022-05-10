function [Pm,search_area]=Bayesian_DOA_root(Y,search_area,etc)
[M,T]=size(Y);
K_hat=length(search_area);

reslu=search_area(2)-search_area(1);
search_mid_left=search_area-reslu/2;
search_mid_right=search_area+reslu/2;


a_search=search_area*pi/180.;
A_theta=exp(-1i*pi*(0:M-1)'*sin(a_search));
%%%%%%%%%%%%%%
a=0.0001;b=0.0001;d=0.01;
maxiter=500;
tol=1e-4;
%initialization
beta=1;
delta=ones(K_hat,1);
%%%%%%%%%%%%%%%

converged = false;
iter = 0;

while ~converged
   iter = iter + 1;
   delta_last = delta;
   
   %calculate mu and Sigma
   Phi=A_theta;
   V_temp= 1/beta*eye(M) +  Phi *diag(delta) * Phi';
   Vinv=inv(V_temp);
   Sigma = diag(delta) -diag(delta) * Phi' * Vinv * Phi *  diag(delta);
   mu = beta * Sigma * Phi' * Y;
   gamma1 = 1 - real(diag(Sigma)) ./ (delta); % delta will change, so move it ahead from updating beta   
      
   %update delta
   temp=sum( mu.*conj(mu), 2) + T*real(diag(Sigma));     
   delta= ( -T+ sqrt(  T^2 + 4*d* real(temp) ) ) / (  2*d   );  
     
   %update beta
   resid=Y-Phi*mu;
   beta=( T*M + (a-1))/( b +  norm(resid, 'fro')^2  +  T / beta * sum(gamma1)     );
 
    % stopping criteria
      erro=norm(delta - delta_last)/norm(delta_last);
    if erro < tol || iter >= maxiter
        converged = true;
    end
    
    
    %%%%%%%%%%%%root-refining
         f=sqrt(sum(mu.*conj(mu),2));
         [~,sort_ind]=sort(f);
         index_amp=sort_ind(end:-1:end-etc+1);

        for j=1:length(index_amp)
            ii=index_amp(j);
            mut=mu(ii,:);
            Sigmat=Sigma(:,ii);
            phi=mut*mut' +  T* Sigmat(ii);
            tempind=[1:K_hat]; tempind(ii)=[];
            Yti=Y- Phi(:, tempind)*mu(tempind,:);
            varphi=  T* Phi(:,tempind)*  (Sigmat(tempind))  -Yti*(mut');
            z1=[1:M-1]';
            c=zeros(M,1);
            c(1)=M*(M-1)/2*phi;
            c(2:end)=z1.*varphi(2:end);
            
           %%%%%%%%% root method 
           ro=roots(c);
           abs_root=abs(ro);
           [~,indmin]=min(abs(abs_root-1));
           angle_cand=  asin(-angle(ro(indmin))/pi)   /pi*180;      
           
           if angle_cand<=search_mid_right(ii) && angle_cand>=search_mid_left(ii)
              search_area(ii)=angle_cand;
              A_theta(:,ii)= exp(-1i*pi*(0:M-1)'*sin(angle_cand*pi/180));               
           end
       
        end    


    
end
Pm=mean(mu.*conj(mu),2);



