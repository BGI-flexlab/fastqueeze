
//text2
// Created by cfy on 2017/5/22.
//instructions:输入的是一条比对后的reads的结果，格式为「chromosome，Pos，Cigar，Strand_Mark」。输出为拆分的数据流。

#include <iostream>
#include <map>
#include <vector>

using namespace std;

struct result_reads
{
    string Chr_name;
    int Pos;
    string Cigar;
    int Strand_Mark;
};

struct aligent   //拆分后的数据流
{
    string Chr_name_pos[2];//染色体+状态（开始F，中间M，结尾E）
    int Pos_relative;  //相对位置，同一条染色体相对上一个reads的尾部的相对位置
    int Pieces;   //记录End的数目
    int Cigar_Mark; //是否有变异，0和1表示
    int Cigar_Num;  //记录同一条染色体中出现的第几个cigar（不带M的那种情况）
    int Position;  //记录位点
    string Cigar_Variation;  //变异情况的编码
    int Strand_Mark_length;   //这个是比对的正反向，0和1表示
};


void readsIn(struct result_reads &result_reads1)
{
    cin>>result_reads1.Chr_name>>result_reads1.Pos>>result_reads1.Cigar>>result_reads1.Strand_Mark;
}

void readsAligent(int &number,vector<aligent>& vec,struct result_reads &result_reads1,struct aligent &AL)
{
    cout<<"输入："<<endl;
    readsIn(result_reads1);
    //判断该Cigar是否存在变异
    bool flag=false;
    string r[5]={"A","T","C","G","N"};
    for(int tt=0;tt<5;tt++)
    {
      if(result_reads1.Cigar.find(r[tt])!=string::npos)
      {flag=true;}
    }

    /************* part 1 *****************/
    //ChrName translates into "ChrName with First and last Pos"
    AL.Chr_name_pos[0]=result_reads1.Chr_name;
    AL.Pieces=1;
    //F表示第一段，M表示中间部分的reads，E表示最后的reads (同一染色体内)
    if(number==0)
    {  //第一条染色体
        AL.Chr_name_pos[1]='F';
    }
    else if(number>0 && (vec[number-1].Chr_name_pos[0]==result_reads1.Chr_name)&&(result_reads1.Pos!=vec[number-1].Position))
    {  //同一条染色体，不同位点
        AL.Chr_name_pos[1]='M';
    }
    else if(number>0 && (vec[number-1].Chr_name_pos[0]==result_reads1.Chr_name) && (result_reads1.Pos == vec[number-1].Position))
    {  //同一条染色体，相同位点
        AL.Chr_name_pos[1]=vec[number-1].Chr_name_pos[1];
        vec[number-1].Pieces++;
        AL.Pieces=vec[number-1].Pieces;
    }
    else if(number>0 &&(vec[number-1].Chr_name_pos[0]!=result_reads1.Chr_name))
    {  //不同染色体
        vec[number-1].Chr_name_pos[1]='E';
        AL.Chr_name_pos[1]='F';
    }


    /************** part 2 *****************/
    //Pos translates into relative position.
    //question：查找上一个结果，如果是同一条染色体的话，写出和上一条染色体之间的相对位置,如果不是同一条染色体的话，写出位置0的相对位置信息。

    vector<int>match;    //保存cigar的分割值，数字
    vector<char>var1;    //保存转换后的突变信息，字符
    bool needPushInt=false;
    int match_num=0;
    AL.Cigar_Num=0;

    for(int temp=0;temp<result_reads1.Cigar.length();temp++)
    {
        char Ci=result_reads1.Cigar[temp];
        if(Ci>='a'&& Ci<='z')
        {
            var1.push_back(Ci-'a'+'A');
            if(needPushInt) {
                match.push_back(match_num);
                match_num = 0;
                needPushInt = false;
            }
        }
        else if((Ci>='A')&&(Ci<='Z'))
        {
            var1.push_back(Ci);
            if(needPushInt)
            {
                match.push_back(match_num);
                match_num=0;
                needPushInt= false;
            }
        }
        else if((Ci>='0') && (Ci<='9'))
        {
            match_num=match_num*10+(Ci-'0');
            needPushInt=true;
        }
    }
    if(needPushInt)
    {
        match.push_back(match_num);
        match_num=0;
        needPushInt= false;
    }
    AL.Cigar_Num=(int)var1.size()/2;
/*
    AL.reads_length=AL.Cigar_Num;   //变异的长度
    for(int N=0;N<match.size();N++)    //计算没有变异的长度
    {
        AL.reads_length=AL.reads_length+match[N];
    }
    cout<<"length: "<<AL.reads_length<<endl;
 */
    AL.Position =result_reads1.Pos;
    if((number==0)||(AL.Chr_name_pos[0]!=vec[number-1].Chr_name_pos[0]))
    {   //染色体的第一条记录，直接取reads的pos
        AL.Pos_relative=result_reads1.Pos;
    }
    else if((number>0)&& (AL.Chr_name_pos[0]==vec[number-1].Chr_name_pos[0]))
    {//同一条染色体的reads比对结果相对上一次的比对结果的位置（在上次的后面Cigar_Relative_Pos位）
        AL.Pos_relative=result_reads1.Pos-(vec[number-1].Position);
    }

    /************* part 3 ********************/
    //Cigars divides into Cigar Mark's Run-Length,Cigars'number,
    // Cigar group's position and Cigar group's variation
    //cigar 的输入格式： 15AT20GC10    其中数字表示match的个数，字符表示突变结果
    if(flag)
    {
        AL.Cigar_Mark=1;  //存在变异
    }
    else AL.Cigar_Mark=0;   //不存在变异

    //变异位置 variations

    map<char,map<char,string> >map2D;
    map<char ,string> map1D;

    map1D['C']="00"; map1D['G']="01"; map1D['T']="10"; map1D['N']="11"; map2D['A']=map1D;

    map1D['G']="00"; map1D['T']="01"; map1D['N']="10"; map1D['A']="11"; map1D.erase('C'); map2D['C']=map1D;

    map1D['T']="00"; map1D['N']="01"; map1D['A']="10"; map1D['C']="11"; map1D.erase('G'); map2D['G']=map1D;

    map1D['N']="00"; map1D['A']="01"; map1D['C']="10"; map1D['G']="11"; map1D.erase('T'); map2D['T']=map1D;

    map1D['A']="00"; map1D['C']="01"; map1D['G']="10"; map1D['T']="11"; map1D.erase('N'); map2D['N']=map1D;
    //varition
    if(flag)  //存在变异
    {
        for(int i=0;i<(var1.size()*2+1);)
        {
            AL.Cigar_Variation.append(map2D[var1[i]][var1[i+1]]);
            i=i+2;
        }
    }
    else  //不存在变异，则无变异记录
    {
        AL.Cigar_Variation="";
    }

    /*********** part 4 **********************/
    //Strand Mark translates into Run-Length
    //匹配的方向，0和1表示
    AL.Strand_Mark_length = result_reads1.Strand_Mark;
    vec.push_back(AL);
    number++;
}

int main()
{
    int number=0;
    int t;
    struct result_reads result_reads1;
    vector<aligent>vec;
    struct aligent AL;

    int kkk=0;
    cout<<"请输出你要测试的数据个数： ";
    cin>>kkk;
    for(int kk=0;kk<kkk;kk++)
    {
        readsAligent(number,vec,result_reads1,AL);
    }

    //测试
    cout<<"******* test vec **********"<<endl;
    for(int k=0;k<number;k++)
    {
        cout<<"Chr_name_pos[0]:       "<<vec[k].Chr_name_pos[0]<<endl;
        cout<<"Chr_name_pos[1]:       "<<vec[k].Chr_name_pos[1]<<endl;
        cout<<"Pieces:                "<<vec[k].Pieces<<endl;
        cout<<"Position:              "<<vec[k].Position<<endl;
        cout<<"Pos_relative:          "<<vec[k].Pos_relative<<endl;
        cout<<"Cigar_Mark:            "<<vec[k].Cigar_Mark<<endl;
        cout<<"Cigar_Num:             "<<vec[k].Cigar_Num<<endl;
        cout<<"Cigar_Variation:       "<<vec[k].Cigar_Variation<<endl;
        cout<<"Strand_Mark_length:    "<<vec[k].Strand_Mark_length<<endl;
        cout<<"******* next reads ********"<<endl;
    }

    //释放vec:
    //vec.erase(vec.begin());
    return 0;
}

/*
sample:
chr1 2 3AT14CG38 1
chr1 3 4AG3CT10 1
chr1 7 10AT20CG29TC 1
chr1 18 10AG39GC17TA2 1
chr1 19 AT18CA27GT28 1
chr1 38 2TC19AG31TC20 1

 * */