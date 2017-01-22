{

    |\      /|  /-------\  |-----\    |-------|
    | \    / |  /       \  |      \   |
    |  \  /  |  /       \  |       \  |
    |   \/   |  /-------\  |       |  |----|
    |        |  /       \  |       /  |
    |        |  /       \  |      /   |
    |        |  /       \  |-----/    |-------|

                ---------  | /----\
                    |      |/      \
                    |      |       |
                    |      |       |
                    |      |       |
                    |      |       |
                ---------  |       |

     ---------  /-------\  /-------\  \-------/
         |      |          |           \-----/
         |      |          |            \---/
         |      \-------\  \-------\     \-/
         |              |          |      |
         |              |          |     |-|
     ---------  \-------/  \-------/     |-|
     
     /-------\  /-------\     /|      /-------\
     |       |  |       |    / |      |      
             |  |       |   /  |      |     
     /-------/  |       |  /   |      /-------\
     |          |       |      |      |       |   
     |          |       |      |      |       |   
     \-------/  \-------/  \-------/  \-------/  
      
}

(*
    Stable version from 17.01.2017
    -Fixed some bugs, including famous bag with Determinant
    -Added some extra types for typedefs
    -Rethinking "arr" and "values" as a property
    -Now powering up function is a CLASS function
*)

unit matrixunit_dev;

interface

const 
  error_msg = 'Error during working with matrix'; //message on error exception
  error_det = 'Determinant of such matrix is unreachable'; //message of unreachable determinant

type
  int = shortint;                          //-128..127  //typedef for int
        //smallint;                          //-32768..32767
        //integer;                           //-2147483648..2147483647
        //longint;                           //-2147483648..2147483647
        //int64;                             //-9223372036854775808..9223372036854775807
        //BigInteger;                        //unlimited
  uint = byte;                             //0..255  //typedef for uint
         //word;                             //0..65535
         //longword;                         //0..4294967295
         //cardinal;                         //0..4294967295
         //uint64                            //0..18446744073709551615
         //BigInteger;                       //unlimited
  float = real;                            //-1.8∙10^308 .. 1.8∙10^308  //typedef for float
          //double;                          //-1.8∙10^308 .. 1.8∙10^308
          //single;                          //-3.4∙10^38 .. 3.4∙10^38
          //decimal                          //-79228162514264337593543950335..79228162514264337593543950335
  bool = boolean;                          //false..true  //typedef for logical bool
  table = aRRay of aRRay of float; //typedef for matrix tables
  line = aRRay of float; //typedef for linear array
  TMatrix = class //main class of TMatrix (matrix)
  private
    farr:table; //array of matrix
    fvalues:line; //array of free coefficients
    class function pow(a, n:int):int; //Resursive function of powering "a" to "n" exponent
    class procedure swap<T>(var a, b:T);//procedure of swapping two variables with ANY type
    //function Formed(a:table; y:int; b:line):table;//+
    //class function Minor(a:table; x, y:int):table; //function of calculating the minor of the matrix
  public
    ///Matrix with float numbers
    property arr:table read farr write farr; 
    ///Linear array with free coefficients
    property values:line read fvalues write fvalues; 
    constructor Create(x, y:uint; flag:bool); //Construstor which allocates memory for "x" lines and "y" columns
    procedure SetLenTable(x, y:uint); //procedure which sets height and width for matrix
    procedure SetLenRates(x:uint); //procedure which sets lenght of linear array with free coefs
    procedure FillByKeyboard(flag:bool); //procedure of filling matrix by keyboard
    procedure FillByFile(fin:TextFile; flag:bool); //procedure of filling matrix by textfile
    procedure FillByMatrix(t:table; l:line); //procedure of filling matrix from other matrix (table and line of free coefs)
    procedure FillRatesByKeyboard; //procedure of filling linear array with free coefs by keyboard
    procedure ShowUp(flag:bool); //procedure of showing matrix content up to screen
    procedure SaveToFile(var fout:TextFile; flag:bool); //procedure of saving matrix content to file
    procedure SaveToMatrix(var t:table; var l:line); //procedure of saving matrix content to table and linear array
    procedure Transpose; //procedure of transposing matrix table
    procedure Rotation90Degrees(flag:bool); //procedure of rotation the matrix for 90 degrees clockwise
    //function Determinant(t:table; n:uint):float; //function of calculating determinant of the matrix table
    //function KramerSolve:list;
    //function GaussSolve:list;}
    //{procedure ShowRatesUp;//+
  end;
  
implementation

const k = -1;

///Recursive function of powering "a" to "n" exponent
class function TMatrix.pow(a, n:int):int; 
begin
  if (n = 1) then //if n equals 1 (exit condition)
    Result:= a //exiting away
  else Result:= a * pow(a, n - 1); //else calling the same function but with other params
end;

///Procedure of swapping two variables with ANY type
class procedure TMatrix.swap<T>(var a, b:T);
begin
  var buf:= a; //saving data from "a"
  a:= b; //replacing data of "a" with data of "b"
  b:= buf; //restoring "buf" data to "b" data
end;

///Construstor which allocates memory for "x" lines and "y" columns.
///Send "true" if you want to fill it by zeros and "false" if you don't need so.
constructor TMatrix.Create(x, y:uint; flag:bool);
var i, j:uint;
begin
  inherited Create; //calling inherited constructor (for object)
  SetLength(farr, x); //setting length of matrix
  for i:= low(farr) to high(farr) do //cycling on matrix
  begin
    SetLength(farr[i], y); //setting length of i-line of matrix
    if (flag) then //if we need filling by zeros
      for j:= low(farr[i]) to high(farr[i]) do //cycling on current line
        farr[i, j]:= 0; //filling by zeros
  end;
  SetLength(fvalues, x); //allocating memory for linear array with free coefs
  if (flag) then //if we need filling by zeros
    for i:= low(fvalues) to high(fvalues) do //cycling on linear array with free coefs
      fvalues[i]:= 0; //filling by zeros
end;

///Procedure which sets height "x" and width "y" for matrix table
procedure TMatrix.SetLenTable(x, y:uint);
const err = 'Cannot allocate memory for matrix table'; //exception name
var i:uint;
begin  
  try
    SetLength(farr, x); //setting height
    for i:= low(farr) to high(farr) do //cycling on matrix table
      SetLength(farr[i], y); //setting width of the current line
  except
    raise Exception.Create(err); //raising following exception
  end;
end;

///Procedure which sets length "x" for linear array with free coefs 
procedure TMatrix.SetLenRates(x:uint);
const err = 'Cannot allocate memory for linear array with free coefs'; //exception name 
begin
  try
    SetLength(fvalues, x); //allcating memory for linear array with free coefs  
  except
    raise Exception.Create(err); //raising following exception
  end;  
end;

///Procedure of filling matrix table by keyboard (STD_INPUT).
///If you want to fill linear array with free coefs too, send "true", else send "false".
procedure TMatrix.FillByKeyboard(flag:bool);
const err = 'Exception during matrix filling'; //exception name
var i, j:uint;
begin
  try
    for i:= low(farr) to high(farr) do //cycling between lines of matrix table
    begin
      for j:= low(farr[i]) to high(farr[i]) do //cycling on current line
        read(farr[i, j]); //reading current item of matrix table
      if (flag) then //if we need to read coefs
        read(fvalues[i]); //reading current item of linear array with free coefs
    end;
  except
    raise Exception.Create(err); //raising following exception
  end;
end;

///Procedure of filling matrix table by file "fin".
///If you want to fill linear array with free coefs too, send "true", else send "false".
procedure TMatrix.FillByFile(fin:TextFile; flag:bool);
const err = 'Unexpected end of file'; //error message
var i, j:uint;
begin
  try
    reset(fin); //resetting file "fin"
    for i:= low(farr) to high(farr) do //cycling of matrix table lines
    begin
      for j:= low(farr[i]) to high(farr[i]) do //cycling of current line
        read(fin, farr[i, j]);
      if (flag) then //if we need to read coefs
        readln(fin, fvalues[i]); //reading current coef
    end;
  except
    raise Exception.Create(err); //raising following exception
  end;
end;

///Procedure of filling matrix by parts (table "t" and linear array with free coefs "l")
///of other matrix. Sending "l" is required.
procedure TMatrix.FillByMatrix(t:table; l:line);
const err = 'Error of filling the matrix or allocating memory'; //error message
var i, j:uint;
begin
  try
    SetLength(farr, Length(t));//allocating memory for table lines
    for i:= low(t) to high(t) do //cycling on table lines
    begin
      SetLength(farr[i], Length(t[i])); //allocating memory for current line
      for j:= low(t[i]) to high(t[i]) do //cycling on current line
        farr[i, j]:= t[i, j]; //copying data to new place
    end;
    SetLength(fvalues, Length(l)); //allocating memory for linear array with free coefs
    for i:= low(l) to high(l) do //cycling on linear array with free coefs
      fvalues[i]:= l[i]; //copying data to new place
  except
    raise Exception.Create(err); //raising following exception
  end;
end;

///Procedure of filling linear array with free coefs by keyboard (STD_INPUT).
procedure TMatrix.FillRatesByKeyboard;
const err = 'Error on filling linear array with free coefs'; //error message
var i:uint;
begin
  try
    for i:= low(fvalues) to high(fvalues) do //cycling on linear array wit free coefs
      read(fvalues[i]); //reading current value
  except
    raise Exception.Create(err); //raising following exception
  end;
end;

///Procedure of showing matrix table content up to screen or console (STD_OUTPUT).
///If you need to show linear array with free coefs too send "true" else send "false".
procedure TMatrix.ShowUp(flag:bool);
const err = 'Error of writing content of matrix to console (STD_OUTPUT)'; //error message
var i, j:uint;
begin
  try
    for i:= low(farr) to high(farr) do //cycling on matrix table lines
    begin
      for j:= low(farr[i]) to high(farr[i]) do //cycling on current line
        write(farr[i, j]:4); //writing out current value
      if (flag) then //if we need to wite out free coefs
        write(fvalues[i]:4); //writing out current free coef-ts
      writeln; //writing over current line in console
    end;
  except
    raise Exception.Create(err); //raising following exception
  end;
end;

///Procedure of saving matrix table content to file "fout".
///If you need to save linear array with free coefs too send "true" else send "false".
procedure TMatrix.SaveToFile(var fout:TextFile; flag:bool);
const err = 'Error of saving matrix content to file'; //error message
var i, j:uint;
begin
  try
    rewrite(fout); //rewriting and erasing output file
    for i:= low(farr) to high(farr) do //cycling on matrix table lines
    begin
      for j:= low(farr[i]) to high(farr[i]) do //cycling on current line
        write(fout, farr[i, j]:4); //saving current value
      if (flag) then //if we need to save free coef-ts
        write(fout, fvalues[i]:4); //saving current free coef
      writeln(fout); //writing over current line in file
    end;
  except
    raise Exception.Create(err); //raising following exception
  end;
end;

///Procedure of saving matrix content to table "t" and linear array "l"
procedure TMatrix.SaveToMatrix(var t:table; var l:line);
const err = 'Error of copying matrix to output table and linear array'; //error message
var i, j:uint;
begin
  try
    SetLength(t, Length(farr)); //allocating memory for output table
    for i:= low(farr) to high(farr) do //cycling on the table
    begin
      SetLength(t[i], Length(farr[i])); //allocating memory for current table line
      for j:= low(farr[i]) to high(farr[i]) do //cycling on current table line
        t[i, j]:= farr[i, j]; //copying current data
    end;
    SetLength(l, Length(fvalues)); //allocating memory for output linear array
    for i:= low(fvalues) to high(fvalues) do //cycling on linear array
      l[i]:= fvalues[i]; //copying current data
  except
    raise Exception.Create(err); //raising following exception
  end;
end;

///Procedure of transposing matrix table
procedure TMatrix.Transpose;
var i, j:uint;
begin
  for i:= low(farr) to high(farr) do //cycling on the matrix table lines
    for j:= i + 1 to high(farr[i]) do //cycling on current line
      swap(farr[i, j], farr[j, i]); //swapping
end;

///Procedure of rotating matrix table for 90 degrees
///Send "true" for rotating clockwise or "false" for rotating anticlockwise
procedure TMatrix.Rotation90Degrees(flag:bool);
var i, j:uint;
begin
  Transpose; //transposing matrix table
  if (flag) then //if we need rotaing clockwise
    for i:= low(farr) to high(farr) do //cycling on the matrix table
      for j:= low(farr[i]) to (Length(farr[i]) div 2) - 1 do //cycling on current line
        swap(farr[i, j], farr[i, high(farr[i]) - j])
  else //else if we need rotating anticlockwise
    for i:= low(farr) to (Length(farr) div 2) - 1 do //cycling on the matrix table
      for j:= low(farr[i]) to high(farr[i]) do //cycling on current line
        swap(farr[i, j], farr[high(farr) - i, j]); //swapping
end;

{///Recursive function of calculating determinant of the table "a".
///Send 0 for "n" as default value (to llok at ALL table).
function TMatrix.Determinant(t:table; n:uint):float;
const err = 'Impossible to calculate the determinant'; //error message
      eps = 0.0000000000001; //machine zero
var det, max, tmp:float;
    i, j, index:uint;
    buf:line;
    b:table;
begin
  if (Length(t) <> Length(t[low(t)])) then //checking equality of height and width of table
    raise Exception.Create(err); //raising following exception
  if (i = 2) then //if we have 2x2 table
  begin
    Result:= t[low(t), low(t[low(t)])] * t[low(t) + 1, low(t[low(t)]) + 1]
             - t[low(t) + 1, low(t[low(t)])] * t[low(t), low(t[low(t)]) + 1]; //calculating determinant simply
    exit; //exiting away
  end;
  det:= 1; //setting determinant as "one"
  index:= low(a); //setting index of a max value item as low index
  for i:= low(t) to high(t) do //cycling on table lines
  if (abs(t[i, low(t[i])]) > abs(t[index, low(t[index])]) then //if we found new abs-value max
    index:= i; //overwriting the index of max
  if (abs(t[index, low(t[index])]) <= eps) then //if maximum equals "zero"
  begin
    Result:= 0; //returning "zero"
    exit; //exiting away
  end;
  if (index <> low(t)) then //if index isn't low index of table lines
  begin
    swap(t[index], t[low(t)]); //swapping values
    det:= -det; //refreshing current determinant
  end;
  max:= t[low(t), low(t[low(t)])]; //getting max abs-value for a special variable
  for j:= low(t[low(t)]) to high(t[low(t)]) do //cycling on columns of low-indexed line of the matrix
  begin
    t[low(t), j]:= t[low(t), j] / max; //dividing current value to "max"
    if (j > low(a[low(a)])) then
    begin
      tmp:= a[low(a), j];
      for i:= low(a) to high(a) do
        a[i, j]:= a[i, j] - tmp * a[i, low(a[i])];
    end;
  end;
  det:= det * max;
  SetLength(b, Length(a) - 1);
  for i:= low(b) to high(b) do
  begin
    SetLength(b[i], Length(a[i]) - 1);
    for j:= low(b[i]) to high(b[i]) do
    begin
      b[i, j]:= a[i + 1, j + 1];
    end;
  end;
  Result:= det * Determinant(b);
end;}

end.

{class function TMatrix.Minor(a:table; x, y:int):table;
var i, j, m:Byte;
    b:table;
    dx, dy:0..1;
begin
  i:= Length(a);
  m:= Length(a[0]);
  SetLength(b, i - 1);
  for i:= low(a) to high(a) do
  begin
    SetLength(mas[i], m - 1);
    for j:= low(a[i]) to high(a[i]) do
    begin
      if (i >= x) then
        dx:= 1
      else
        dx:= 0;
      if (j >= y) then
        dy:= 1
      else
        dy:= 0;
      b[i, j]:= a[i + dx, j + dy];
    end;
  end;
  for i:= low(b) to high(b) do
  begin
    for j:= low(b[i]) to high(b[i]) do
      write(b[i, j]:4);
    writeln;
  end;
  Result:= b;
end;}

{function TMatrix.Formed(a:table; y:int; b:list):table;
var i:int;
begin
  if (Length(a) <> Length(b)) then
    raise Exception.Create('Размер исходного столбца не совпадает с размером заменителя');
  for i:= low(a) to high(a) do
    a[i, y]:= b[i];
  Result:= a;
end;}

{procedure TMatrix.ShowRatesUp;
var i:unsigned;
begin
  try
    for i:= low(value) to high(value) do
      write(value[i]:4);
  except
    raise Exception.Create('Ошибка при выводе коэф-тов на экран');
  end;
end;}