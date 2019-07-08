classdef Operator<handle
%- Author:石凯元
%- Time: 08 Jul 2019
%- Please follow GPL License using the source code
    properties
        indexs;operator;
    end
    methods
        function self = Operator(center,indexs,operator)
            if strcmp(class(operator),'Operator')
                self=operator;
                return
            end
            if length(center)==0
                num=max(sum(operator~=0,2));
                len=length(operator);
                self.indexs=cell(num,1);
                self.operator=cell(num,1);
                self.indexs(:)={ones(len,1)};
                self.operator(:)={zeros(len,1)};
                [row,col]=find(operator);
                for i=1:length(row)
                    curRow=row(i);
                    curCol=col(i);
                    cur=1;
                    while self.operator{cur}(curRow)~=0
                        cur=cur+1;
                    end
                    self.indexs{cur}(curRow)=curCol;
                    self.operator{cur}(curRow)=operator(curRow,curCol);
                end
            else
                self.indexs=indexs;
                for i=1:length(indexs)
                    self.operator{i}=operator(sub2ind(size(operator),center,indexs{i}));
                end
            end
        end
        function newOp = plus(op1,op2)
            newOp=Operator([],{},[]);
            newOp.indexs=op1.indexs;
            for i=1:length(op1.indexs)
                newOp.operator{i}=op1.operator{i}+op2.operator{i};
            end
        end
        function newOp = minus(op1,op2)
            newOp=Operator([],{},[]);
            newOp.indexs=op1.indexs;
            for i=1:length(op1.indexs)
                newOp.operator{i}=op1.operator{i}-op2.operator{i};
            end
        end
        function newOp = times(array,op)
            newOp=Operator([],{},[]);
            newOp.indexs=op.indexs;
            for i=1:length(op.indexs)
                newOp.operator{i}=array.*op.operator{i};
            end
        end
        function newarray = mtimes(op,array)
            newarray=op.operator{1}.*array(op.indexs{1});
            for i=2:length(op.indexs)
                newarray=newarray+op.operator{i}.*array(op.indexs{i});
            end
        end
        function gpuIn(self)
            self=self;
        end
        function gpuOut(self)
            self=self;
        end
    end
end