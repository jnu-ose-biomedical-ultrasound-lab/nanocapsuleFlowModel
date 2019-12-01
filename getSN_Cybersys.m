%% get the size and number of RBC aggregation
% A = AGG_p;
function [Agg_size Agg_No AG a RBC_2]= getSN_Cybersys2(A)
 B=sort([A(:,1); A(:,2)]);                                         %每一列从小到大排列，找出有多少个RBC
B_1=unique(A);                                                    %找出B的唯一值
C=B(find(diff(B)==0));                                         %找出重复的RBC的位置对应的值
if length(C)==0
    No_2=length(B_1)/2;
    AG=[2*ones(1,No_2)];
    Agg_size=2;
    Agg_No=No_2;
    a={};
    RBC_2=B_1;
  
    else
    sz=size(A,1);
    kk=1;
    To_1=[];
    for ii=1:length(C)                                                 %重复的个数
        inx2=[];
        D=find(A==C(ii));                                            %把A的两列合成一列，找到重复的位置
        inx1=D-sz*((D-sz)>0);                                     %找到在A的位置在A的列长之内的值
        g_e1=unique(A(inx1,:));
        for jj=1:length(g_e1)                                        %
            D2=find(A==g_e1(jj));
            inxh=D2-sz*((D2-sz)>0);
            inx2=[inx2; inxh];
        end
            inx=unique([inx1;inx2]);
            a{ii}=unique(A(inx,:));
            To_1=[To_1;a{ii}];
            if (ii>1 & kk<ii+1)
                for kk=1:ii-1
                    com_1=intersect(a{ii},a{kk});% too slow here
                    if length(com_1)>0
                        a{kk}= unique([a{ii};a{kk}]);
                        a{ii}=[];
                    end
                end
        end
     end
        Csize = cellfun('length', a);
    anyzero = find(Csize ==0);
    Csize(anyzero)=[];
    To_1=unique(To_1);
    RBC_2=setdiff(B_1,To_1);
    No_2=length(RBC_2)/2;
    AG=[Csize 2*ones(1,No_2)];
    Agg_size=mean(AG);
    Agg_No=length(AG);
end
end
