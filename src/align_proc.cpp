//
// Created by cfy on 2017/6/7.
//
#include "encode.h"

using namespace std;


/**********************  encode  ************************************/
void maplist(map<char,map<char,string> >& map2D)
{
    //变异位置 variations
    //map<char,map<char,string> >map2D;
    map<char ,string> map1D;

    map1D['C']="00"; map1D['G']="01"; map1D['T']="10"; map1D['N']="11"; map2D['A']=map1D;

    map1D['G']="00"; map1D['T']="01"; map1D['N']="10"; map1D['A']="11"; map1D.erase('C'); map2D['C']=map1D;

    map1D['T']="00"; map1D['N']="01"; map1D['A']="10"; map1D['C']="11"; map1D.erase('G'); map2D['G']=map1D;

    map1D['N']="00"; map1D['A']="01"; map1D['C']="10"; map1D['G']="11"; map1D.erase('T'); map2D['T']=map1D;

    map1D['A']="00"; map1D['C']="01"; map1D['G']="10"; map1D['T']="11"; map1D.erase('N'); map2D['N']=map1D;
}

void readsIn(struct result_reads &result_reads1)
{
    cin>>result_reads1.Chr_name>>result_reads1.Pos>>result_reads1.Cigar>>result_reads1.Strand_Mark;
}

void align_encode(int &number,vector<aligent>& vec,map<string,map<int,Chr> >& chr)
{
    map<char,map<char,string> >map2D;
    struct result_reads result_reads1;
    struct aligent AL;
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


    //ChrName translates into "ChrName with First and last Pos"
    AL.Chr_name_pos[0]=result_reads1.Chr_name;
    //initialization
    AL.Pieces=1;
    chr[AL.Chr_name_pos[0]][0].Pieces=1;

    //F——first, M——middle, E——end;
    if(number==0)
    {  //第一条染色体
        AL.Chr_name_pos[1]='F';
        chr[AL.Chr_name_pos[0]][0].Mark='F';
        chr[AL.Chr_name_pos[0]][0].Position=result_reads1.Pos;
        chr[AL.Chr_name_pos[0]][0].Pieces=1;
    }
    else if(number>0 && (vec[number-1].Chr_name_pos[0]==result_reads1.Chr_name)&&(result_reads1.Pos!=vec[number-1].Position))
    {  //同一条染色体，不同位点
        //AL.Chr_name_pos[1]='M';
        if((vec[number-1].Chr_name_pos[1]=="F")||(vec[number-1].Chr_name_pos[1]=="M"))
        {
            AL.Chr_name_pos[1]='L';
        }
        else if(vec[number-1].Chr_name_pos[1]=="L")
        {
            vec[number-1].Chr_name_pos[1]='M'; AL.Chr_name_pos[1]='L';
        }
        chr[AL.Chr_name_pos[0]][1].Mark='L';
        chr[AL.Chr_name_pos[0]][1].Pieces=1;
        chr[AL.Chr_name_pos[0]][1].Position=result_reads1.Pos;
    }
    else if(number>0 && (vec[number-1].Chr_name_pos[0]==result_reads1.Chr_name) && (result_reads1.Pos == vec[number-1].Position))
    {  //同一条染色体，相同位点
        //AL.Chr_name_pos[1]=vec[number-1].Chr_name_pos[1];
        AL.Chr_name_pos[1]='L';
        vec[number-1].Pieces++;
        AL.Pieces=vec[number-1].Pieces;
        chr[AL.Chr_name_pos[0]][1].Position=result_reads1.Pos;
        chr[AL.Chr_name_pos[0]][1].Mark='L';
        chr[AL.Chr_name_pos[0]][1].Pieces=vec[number-1].Pieces;

        if(vec[number-1].Chr_name_pos[1]=="F")  //F 的 piece
        {
            chr[AL.Chr_name_pos[0]][0].Pieces=vec[number-1].Pieces;
            chr[AL.Chr_name_pos[0]][1].Pieces=vec[number-1].Pieces;  //同样的reads出现在头尾
        }
    }
    else if(number>0 &&(vec[number-1].Chr_name_pos[0]!=result_reads1.Chr_name))
    {  //不同染色体
        //vec[number-1].Chr_name_pos[1]='E';
        vec[number-1].Chr_name_pos[1]='L';
        chr[vec[number-1].Chr_name_pos[0]][1].Pieces=vec[number-1].Pieces;  //End 的 piece

        AL.Chr_name_pos[1]='F';
        chr[vec[number-1].Chr_name_pos[0]][1].Mark='L';  //
        chr[vec[number-1].Chr_name_pos[0]][1].Position=vec[number-1].Position;
        chr[AL.Chr_name_pos[0]][0].Mark='F';
        chr[AL.Chr_name_pos[0]][0].Position=result_reads1.Pos;
    }

    //Pos translates into relative position.
    //question：查找上一个结果，如果是同一条染色体的话，写出和上一条染色体之间的相对位置,如果不是同一条染色体的话，写出位置0的相对位置信息。

    //vector<int>match;    //保存cigar的分割值，数字
    vector<char>var1;    //保存转换后的突变信息，字符
    //bool needPushInt=false;
    int match_num=0;
    AL.Cigar_Num=0;

    for(int temp=0;temp<result_reads1.Cigar.length();temp++)
    {
        char Ci=result_reads1.Cigar[temp];
        if(Ci>='a'&& Ci<='z')
        {
            var1.push_back(Ci-'a'+'A');
        }
        else if((Ci>='A')&&(Ci<='Z'))
        {
            var1.push_back(Ci);
        }
    }

    AL.Cigar_Num=(int)var1.size()/2;

    match_num=0;
    for(int i=0;i<result_reads1.Cigar.length();i++)
    {
        if(result_reads1.Cigar[i]>='a'&&result_reads1.Cigar[i]<='z')
        {
            AL.Cigra_relative_Pos.push_back(0);
            i++;  //出现一个字符接着跳过下一个字符；
        }
        else if(result_reads1.Cigar[i]>='A'&&result_reads1.Cigar[i]<='Z')
        {
            AL.Cigra_relative_Pos.push_back(0);
            i++;  //出现一个字符接着跳过下一个字符；
        }
        else if(result_reads1.Cigar[i]>='0'&&result_reads1.Cigar[i]<='9')
        {
            match_num=match_num*10+(result_reads1.Cigar[i]-'0');
            if((result_reads1.Cigar[i+1]>='A'&& result_reads1.Cigar[i+1]<='Z')||(i==result_reads1.Cigar.length()-1))
            { AL.Cigra_relative_Pos.push_back(match_num); match_num=0;}
        }
    }

    /*
    AL.reads_length=AL.Cigar_Num;   //变异的长度
    for(int N=0;N<AL.match.size();N++)    //计算没有变异的长度
    { AL.reads_length=AL.reads_length+AL.match[N];
    cout<<AL.match[N]<<endl;}
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

    //Cigars divides into Cigar Mark's Run-Length,Cigars'number,
    // Cigar group's position and Cigar group's variation
    //cigar 的输入格式： 15AT20GC10    其中数字表示match的个数，字符表示突变结果
    if(flag)
    {
        AL.Cigar_Mark=1;  //存在变异
    }
    else AL.Cigar_Mark=0;   //不存在变异

    //变异位置 variations

    maplist(map2D);
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


    //Strand Mark translates into Run-Length
    //匹配的方向，0和1表示
    AL.Strand_Mark_length = result_reads1.Strand_Mark;
    vec.push_back(AL);
    number++;
}

void vec_erase(vector<aligent>& vec)
{
    vec.erase(vec.begin());
}




/*************  decode   ************************/

void Reference1(map<string,string>& reference)
{
    reference["chr1"]="AATTAAAAAAATCGATCGATCCTAGATCGGACAACTCGATCGATCGGATCATCGATCGATCGATCCTAGATCGGACAACTCGATCGATCGGATCATCGATCGATCGATCCTAGATCGGACAACTCGATCGATCGGATCATCGATCGATCGATCCTAGATCGGACAACTCGATCGATCGG";
    reference["chr2"]="ATCATCGATCGATCGATCGATCGATCGGATCATCGATCGATCGATCCTAGATCGGACAACTCGATCGATCGGATCATCGATCGATCGATCCTAGATCGGACAACTCGATCGATCGGATCATCGATCGATCGATCCTAGATCGGACAACTCGATCGATCGGATCATCGATCGATCGATCCTAGATCGGACAACTCGATCGATCGG";
    reference["chr3"]="ATCATCGATCGATCGATCGATCGATCGGATCATCGATCGATCGATCCTAGATCGGACAACTCGATCGATCGGATCATCGATCGATCGATCCTAGATCGGACAACTCGATCGATCGGATCATCGATCGATCGATCCTAGATCGGACAACTCGATCGATCGG";
    reference["chr4"]="ATCATCGATCGATCGATCATCGATCGATCGATCCTAGATCGGACAACTATCATCGATCGATCGATCCATCATCGATCGATCGATCCTAGATCGGACAACTCGATCGATCGGTAGATCGGACAACTCGATCGATCGGCGATCGATCGGATCGATCGATCGG";
}
void align_decode(vector<aligent>& vec_decode,map<string,map<int,Chr> >& chr)
{
    aligent AL_decode;   //中间量
    result_reads result_reads2;
    map<char,map<char,string> >map2D;
    map<string,string>reference;

    //initialization
    AL_decode.Chr_name_pos[0]=vec_decode[0].Chr_name_pos[0];
    AL_decode.Chr_name_pos[1]=vec_decode[0].Chr_name_pos[1];
    AL_decode.Pieces=vec_decode[0].Pieces;
    AL_decode.Position=vec_decode[0].Position;
    AL_decode.Pos_relative=vec_decode[0].Pos_relative;
    AL_decode.Cigar_Mark=vec_decode[0].Cigar_Mark;
    AL_decode.Cigar_Num=vec_decode[0].Cigar_Num;
    AL_decode.Cigar_Variation=vec_decode[0].Cigar_Variation;
    AL_decode.Strand_Mark_length=vec_decode[0].Strand_Mark_length;
    AL_decode.Cigra_relative_Pos=vec_decode[0].Cigra_relative_Pos;

    result_reads2.Chr_name=AL_decode.Chr_name_pos[0];
    result_reads2.Pos=AL_decode.Position;
    result_reads2.Strand_Mark=AL_decode.Strand_Mark_length;

    //Cigar的还原
    maplist(map2D);
    Reference1(reference);
    int i=0;
    for(int j=0;j<AL_decode.Cigra_relative_Pos.size();j++)
    {
        if(AL_decode.Cigra_relative_Pos[j]!=0)   //cigar的相对位置不为0的时候是match的序列长度
        {
            AL_decode.Position=AL_decode.Position+AL_decode.Cigra_relative_Pos[j];   //Cigar位置变化
            string s="";
            s+=to_string(AL_decode.Cigra_relative_Pos[j]);   //int转为string
            result_reads2.Cigar.append(s);
        }
        else if(AL_decode.Cigra_relative_Pos[j]==0) {  //Cigar相对位置为0的时候是突变位点，要还原
            switch (AL_decode.Cigar_Variation[i]) {   //变异位点的还原
                case '0': {
                    if (AL_decode.Cigar_Variation[i + 1] == '0') {
                        if (reference[AL_decode.Chr_name_pos[0]][AL_decode.Position] == 'C')
                        { result_reads2.Cigar.append("CG"); }
                        else if (reference[AL_decode.Chr_name_pos[0]][AL_decode.Position] == 'G')
                        { result_reads2.Cigar.append("GT"); }
                        else if (reference[AL_decode.Chr_name_pos[0]][AL_decode.Position] == 'T')
                        { result_reads2.Cigar.append("TN"); }
                        else if (reference[AL_decode.Chr_name_pos[0]][AL_decode.Position] == 'N')
                        { result_reads2.Cigar.append("NA"); }
                        else if (reference[AL_decode.Chr_name_pos[0]][AL_decode.Position] == 'A')
                        { result_reads2.Cigar.append("AC"); }
                    }
                    else if (AL_decode.Cigar_Variation[i + 1] == '1') {
                        if (reference[AL_decode.Chr_name_pos[0]][AL_decode.Position] == 'C')
                        { result_reads2.Cigar.append("CT"); }
                        else if (reference[AL_decode.Chr_name_pos[0]][AL_decode.Position] == 'G')
                        { result_reads2.Cigar.append("GN"); }
                        else if (reference[AL_decode.Chr_name_pos[0]][AL_decode.Position] == 'T')
                        { result_reads2.Cigar.append("TA"); }
                        else if (reference[AL_decode.Chr_name_pos[0]][AL_decode.Position] == 'N')
                        { result_reads2.Cigar.append("NC"); }
                        else if (reference[AL_decode.Chr_name_pos[0]][AL_decode.Position] == 'A')
                        { result_reads2.Cigar.append("AG"); }
                    }
                    break;
                }

                case '1': {
                    if (AL_decode.Cigar_Variation[i + 1] == '0') {
                        if (reference[AL_decode.Chr_name_pos[0]][AL_decode.Position] == 'C')
                        { result_reads2.Cigar.append("CN"); }
                        else if (reference[AL_decode.Chr_name_pos[0]][AL_decode.Position] == 'G')
                        { result_reads2.Cigar.append("GA"); }
                        else if (reference[AL_decode.Chr_name_pos[0]][AL_decode.Position] == 'T')
                        { result_reads2.Cigar.append("TC"); }
                        else if (reference[AL_decode.Chr_name_pos[0]][AL_decode.Position] == 'N')
                        { result_reads2.Cigar.append("NG"); }
                        else if (reference[AL_decode.Chr_name_pos[0]][AL_decode.Position] == 'A')
                        { result_reads2.Cigar.append("AT"); }
                    }
                    else if (AL_decode.Cigar_Variation[i + 1] == '1') {
                        if (reference[AL_decode.Chr_name_pos[0]][AL_decode.Position] == 'C')
                        { result_reads2.Cigar.append("CA"); }
                        else if (reference[AL_decode.Chr_name_pos[0]][AL_decode.Position] == 'G')
                        { result_reads2.Cigar.append("GC"); }
                        else if (reference[AL_decode.Chr_name_pos[0]][AL_decode.Position] == 'T')
                        { result_reads2.Cigar.append("TG"); }
                        else if (reference[AL_decode.Chr_name_pos[0]][AL_decode.Position] == 'N')
                        { result_reads2.Cigar.append("NT"); }
                        else if (reference[AL_decode.Chr_name_pos[0]][AL_decode.Position] == 'A')
                        { result_reads2.Cigar.append("AN"); }
                    }
                    break;
                }
                default: break;
            }
            i=i+2;  //variation每次取两位的
            AL_decode.Position++;
        }
    }
    cout<<"name: "<<result_reads2.Chr_name<<endl;
    cout<<"Pos: "<<result_reads2.Pos<<endl;
    cout<<"Cigar: "<<result_reads2.Cigar<<endl;
}

