---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.11.2
  kernelspec:
    display_name: Julia 1.7.0-DEV
    language: julia
    name: julia-1.7
---

<!-- #region -->
# Expression problem

* https://miguelraz.github.io/blog/dispatch/index.html

以下のコードは

* https://eli.thegreenplace.net/2016/the-expression-problem-and-its-solutions/

より。

__C++__

```C++
class Expr {
public:
  virtual std::string ToString() const = 0;
  virtual double Eval() const = 0;
};
```

```C++
class Constant : public Expr {
public:
  Constant(double value) : val(value) {}

  std::string ToString() const {
    std::ostringstream ss;
    ss << value_;
    return ss.str();
  }

  double Eval() const {
    return value_;
  }

private:
  double value_;
};
```

```C++
class BinaryPlus : public Expr {
public:
  BinaryPlus(const Expr& lhs, const Expr& rhs) : lhs_(lhs), rhs_(rhs) {}

  std::string ToString() const {
    return lhs_.ToString() + " + " + rhs_.ToString();
  }

  double Eval() const {
    return lhs_.Eval() + rhs_.Eval();
  }

private:
  const Expr& lhs_;
  const Expr& rhs_;
};
```

__Haskell__

```Haskell
module Expressions where

data Expr = Constant Double
          | BinaryPlus Expr Expr

stringify :: Expr -> String
stringify (Constant c) = show c
stringify (BinaryPlus lhs rhs) = stringify lhs
                                ++ " + "
                                ++ stringify rhs

evaluate :: Expr -> Double
evaluate (Constant c) = c
evaluate (BinaryPlus lhs rhs) = evaluate lhs + evaluate rhs
```
<!-- #endregion -->

## 値の四則演算を行えるモジュールO

以下をモジュールOで定義する。

* 式の型 `Expression`
  * 定数の型 `Constant`
    * 定数の値を得る函数 `value`
  * 二項演算達の型 `Plus`, `Minus`, `Mul`, `Dic`
    * 二項演算の左辺と右辺を得る函数 `lhs`, `rhs`
* これらを文字列に変換する函数 `stringify`
* そららの値を評価する函数 `evaluate`

```julia
module O

abstract type Expression end

struct Constant{C} <: Expression value::C end
value(c::Constant) = getfield(c, :value)
stringify(c::Constant) = string(value(c))
evaluate(c::Constant) = value(c)

abstract type BinOp <: Expression end
lhs(bo::BinOp) = getfield(bo, :lhs)
rhs(bo::BinOp) = getfield(bo, :rhs)

const binop_list = ((:Plus, :+), (:Minus, :-), (:Mul, :*), (:Div, :/))
for (S, op) in binop_list
    @eval begin
        struct $S{L<:Expression, R<:Expression} <: BinOp lhs::L; rhs::R end
        stringify(bo::$S) = "(" * stringify(lhs(bo)) * " $($op) " * stringify(rhs(bo)) * ")"
        evaluate(bo::$S) = $op(evaluate(lhs(bo)), evaluate(rhs(bo)))
    end
end

# デフォルトでの表示の仕方
Base.show(io::IO, expr::Expression) = print(io, stringify(expr))

end
```

```julia
a = O.Constant(2)
b = O.Constant(3)
c = O.Constant(4)
d = O.Constant(5)
e = O.Constant(6)
a, b, c, d, e
```

```julia
expr1 = O.Div(O.Plus(O.Minus(c, a), O.Mul(b, e)), d)
```

```julia
O.evaluate(expr1)
```

## モジュールOを変更することなく、モジュールPで型と函数を追加

モジュールPでは以下をO.Expression型に追加する。

* 函数呼び出しの型 `FunctionCall`
  * 函数呼び出しの函数部分を取り出す函数 `fun`
  * 函数呼び出しの引数部分を取り出す函数 `arg`
  * 函数呼び出しを文字列に変換するメソッド `O.stringify`
  * 函数呼び出した結果を得るメソッド `O.evaluate`
* 変数の型 `Variable`
  * 変数の名前を得る函数 `name`
  * 変数を文字列に変換するメソッド `O.stringify`
* 変数の値のリストの型 `ValueList`
  * 変数の値のリストの内部実装を取り出す函数 `parent`
  * 変数の値のリストを文字列に変換するメソッド `O.stringify`
  * 変数の値のリストに従って式を評価する函数 `evaluate`
  * 変数の値のリストに従って変数に式を代入した式を作る函数 `substitute`

```julia
module P

### 函数呼び出し機能を追加

using ..O: O, Expression

struct FunctionCall{F, X<:Expression} <: Expression f::F; x::X end
fun(fc::FunctionCall) = getfield(fc, :f)
arg(fc::FunctionCall) = getfield(fc, :x)
O.stringify(fc::FunctionCall) = string(fun(fc)) * "(" * O.stringify(arg(fc)) * ")"
O.evaluate(fc::FunctionCall) = fun(fc)(O.evaluate(arg(fc)))

### 変数を扱う機能を追加

using ..O: O, Expression, Constant, BinOp, lhs, rhs

struct Variable <:Expression name::Symbol end
name(var::Variable) = getfield(var, :name)
O.stringify(var::Variable) = string(name(var))

struct ValueList{T<:NamedTuple} list::T end
ValueList() = ValueList((;))
ValueList(; p...) = ValueList((; p...))
parent(list::ValueList) = getfield(list, :list)
names(list::ValueList) = keys(parent(list))
value(list::ValueList, name::Symbol) = getproperty(parent(list), name)
O.stringify(list::ValueList) = "ValueList" * string(parent(list))

evaluate(var::Variable, list::ValueList) = value(list, name(var))
evaluate(fc::FunctionCall, list::ValueList) = fun(fc)(evaluate(arg(fc), list))
evaluate(c::Constant, list::ValueList) = O.value(c)
for (S, op) in O.binop_list
    @eval function evaluate(bo::O.$S, list::ValueList=ValueList())
        $op(evaluate(lhs(bo), list), evaluate(rhs(bo), list))
    end
end

substitute(var::Variable, list::ValueList) =
    name(var) ∈ names(list) ? Constant(value(list, name(var))) : var
substitute(fc::FunctionCall, list::ValueList) =
    FunctionCall(fun(fc), substitute(arg(fc), list))
substitute(c::Constant, list::ValueList) = c
for (S, op) in O.binop_list
    @eval substitute(bo::O.$S, list::ValueList) =
        O.$S(substitute(lhs(bo), list), substitute(rhs(bo), list))
end

# デフォルトでの表示の仕方
Base.show(io::IO, list::ValueList) = print(io, O.stringify(list))

end
```

```julia
P.ValueList(u=2, v=3, w=4, x=5, y=6)
```

```julia
u = P.Variable(:u)
v = P.Variable(:v)
w = P.Variable(:w)
x = P.Variable(:x)
y = P.Variable(:y)
u, v, w, x, y
```

```julia
sinpi_x_div_y = P.FunctionCall(sinpi, O.Div(x, y))
```

```julia tags=[]
P.evaluate(sinpi_x_div_y, P.ValueList(x = 1, y = 6))
```

```julia tags=[]
P.substitute(sinpi_x_div_y, P.ValueList(x = 1, y = 6))
```

```julia
P.evaluate(expr1, P.ValueList())
```

```julia
expr2 = O.Div(O.Plus(O.Minus(w, u), O.Mul(v, y)), x)
```

```julia
P.evaluate(expr2, P.ValueList(u=2, v=3, w=4, x=5, y=6))
```

```julia
expr3 = P.substitute(expr2, P.ValueList(u=2, v=3, w=4))
```

```julia
expr4 = P.substitute(expr3, P.ValueList(x = 5, y = 6))
```

```julia
O.evaluate(expr4)
```

```julia
P.evaluate(expr3, P.ValueList(x = 5, y = 6))
```

<!-- #region -->
## Visitor patternでモジュールOを全面的に書き直し

以下のコードは

* https://eli.thegreenplace.net/2016/the-expression-problem-and-its-solutions/

より。

__C++__

```C++
class ExprVisitor {
public:
  virtual void VisitConstant(const Constant& c) = 0;
  virtual void VisitBinaryPlus(const BinaryPlus& bp) = 0;
};
```

```C++
class Expr {
public:
  virtual void Accept(ExprVisitor* visitor) const = 0;
};
```

```C++
class Constant : public Expr {
public:
  Constant(double value) : val(value) {}

  void Accept(ExprVisitor* visitor) const {
    visitor->VisitConstant(*this);
  }

  double GetValue() const {
    return value_;
  }

private:
  double value_;
};

// ... similarly, BinaryPlus would have
//
//    void Accept(ExprVisitor* visitor) const {
//      visitor->VisitBinaryPlus(*this);
//    }
//
// ... etc.
```

```C++
class Evaluator : public ExprVisitor {
public:
  double GetValueForExpr(const Expr& e) {
    return value_map_[&e];
  }

  void VisitConstant(const Constant& c) {
    value_map_[&c] = c.GetValue();
  }

  void VisitBinaryPlus(const BinaryPlus& bp) {
    bp.GetLhs().Accept(this);
    bp.GetRhs().Accept(this);
    value_map_[&bp] = value_map_[&(bp.GetLhs())] + value_map_[&(bp.GetRhs())];
  }

private:
  std::map<const Expr*, double> value_map_;
};
```
<!-- #endregion -->

以上のコードを参考にしながら、Visitor patternに従って、モジュールOを完全に作り直してしまいましょう。 

新モジュール名をQとする。

以下のモジュールQはモジュールOよりも圧倒的に複雑になってしまっている。

* 式の型 `Expression`
  * 定数の型 `Constant`
    * 定数の値を得る函数 `value`
  * 二項演算の型 `Plus`, `Minus`, `Mul`, `Div`
    * 二項演算の左辺と右辺を得る函数 `lhs`, `rhs`
* 式訪問者の型 `ExprVisitor`
  * 式訪問者の記憶を取り出す函数 `memory`
  * 式訪問者が式を訪れて仕事をした結果を得る函数 `(::ExprVisitor)(::Expression)`
  * 式訪問者が式を訪れて記憶(思い出)を得る函数 `visit`
  * 式を文字列化する者の型 `Stringifier`
  * 式の評価者の型 `Evaluator`

```julia
module Q

abstract type Expression end

struct Constant{T} <: Expression value::T end
value(x::Constant) = getfield(x, :value)

abstract type BinOp <: Expression end
lhs(bo::BinOp) = getfield(bo, :lhs)
rhs(bo::BinOp) = getfield(bo, :rhs)

const binop_list = ((:Plus, :+), (:Minus, :-), (:Mul, :*), (:Div, :/))
for (S, op) in binop_list
    @eval struct $S{L<:Expression, R<:Expression} <: BinOp lhs::L; rhs::R end
end

abstract type ExprVisitor end
memory(visitor::ExprVisitor) = getfield(visitor, :memory)
function (visitor::ExprVisitor)(expr::Expression) # 多重ディスパッチを利用
    visit(visitor, expr)
    memory(visitor)[expr]
end

struct Evaluator <: ExprVisitor memory::Dict{Expression, Any} end
Evaluator() = Evaluator(Dict{Expression, Any}())

# 多重ディスパッチを利用
visit(visitor::Evaluator, c::Constant) = memory(visitor)[c] = value(c)
for (S, op) in binop_list
    @eval function visit(visitor::Evaluator, bo::$S)
        visit(visitor, lhs(bo))
        visit(visitor, rhs(bo))
        memory(visitor)[bo] = $op(memory(visitor)[lhs(bo)], memory(visitor)[rhs(bo)])
    end
end

struct Stringifier <: ExprVisitor memory::Dict{Expression, String} end
Stringifier() = Stringifier(Dict{Expression, String}())

# 多重ディスパッチを利用
visit(visitor::Stringifier, c::Constant) = memory(visitor)[c] = string(value(c))
for (S, op) in binop_list
    @eval function visit(visitor::Stringifier, bo::$S)
        visit(visitor, lhs(bo))
        visit(visitor, rhs(bo))
        memory(visitor)[bo] =
            "(" * memory(visitor)[lhs(bo)] * " $($op) " * memory(visitor)[rhs(bo)] * ")"
    end
end

# デフォルトでの表示の仕方
Base.show(io::IO, expr::Expression) = print(io, Stringifier()(expr))

end
```

```julia
expr = Q.Plus(Q.Constant(1), Q.Constant(2))
```

```julia
stringifier = Q.Stringifier()
Q.visit(stringifier, expr)
Q.memory(stringifier)
```

```julia
Q.Stringifier()(expr)
```

```julia
evaluator = Q.Evaluator()
Q.visit(evaluator, expr)
Q.memory(evaluator)
```

```julia
Q.Evaluator()(expr)
```

<!-- #region -->
## Visitor patternでの機能の追加

以下のコードは

* https://eli.thegreenplace.net/2016/the-expression-problem-and-its-solutions/

より。

__C++__

```C++
class Evaluator : virtual public ExprVisitor {
  // .. the rest is the same
};
```

```C++
// This is the new ("extended") expression we're adding.
class FunctionCall : public Expr {
public:
  FunctionCall(const std::string& name, const Expr& argument)
      : name_(name), argument_(argument) {}

  void Accept(ExprVisitor* visitor) const {
    ExprVisitorWithFunctionCall* v =
        dynamic_cast<ExprVisitorWithFunctionCall*>(visitor);
    if (v == nullptr) {
      std::cerr << "Fatal: visitor is not ExprVisitorWithFunctionCall\n";
      exit(1);
    }
    v->VisitFunctionCall(*this);
  }

private:
  std::string name_;
  const Expr& argument_;
};
```

```C++
class ExprVisitorWithFunctionCall : virtual public ExprVisitor {
public:
  virtual void VisitFunctionCall(const FunctionCall& fc) = 0;
};
```

```C++
class EvaluatorWithFunctionCall : public ExprVisitorWithFunctionCall,
                                  public Evaluator {
public:
  void VisitFunctionCall(const FunctionCall& fc) {
    std::cout << "Visiting FunctionCall!!\n";
  }
};
```
<!-- #endregion -->

Visitor patternでモジュールOに函数呼び出しの機能を追加してみる。 

多重ディスパッチのおかげで、上のＣ＋＋版と比較すると圧倒的にシンプルになっている。 (函数呼び出しのコードのみを比較せよ)

しかし、モジュールPにおける次の函数呼び出しのコードと比較するとロジックは複雑で込み入ったものになってしまう。

```Julia
struct FunctionCall{F, X<:Expression} <: Expression f::F; x::X end
fun(fx::FunctionCall) = getfield(fx, :f)
arg(fx::FunctionCall) = getfield(fx, :x)
O.stringify(fx::FunctionCall) = string(fun(fx)) * "(" * O.stringify(arg(fx)) * ")"
O.evaluate(fx::FunctionCall) = fun(fx)(O.evaluate(arg(fx)))
```

モジュールRでは以下を追加する。

* 函数呼び出しの型 `FunctionCall`
  * 函数呼び出しの函数部分を取り出す函数 `fun`
  * 函数呼び出しの引数部分を取り出す函数 `arg`
* 式を文字列に変換する者が函数呼び出しを訪れて記憶を得るメソッド `Q.visit`
* 式評価者が函数呼び出しを訪れて記憶を得るメソッド `Q.visit`
* 変数の型 `Variable`
  * 式を文字列に変換する者が変数を訪れて記憶を得るメソッド `Q.visit`
* 変数の値のリストの型 `ValueList`
  * 変数の値のリストを文字列に変換する函数 `stringify`
* 変数の値のリストの下で、式訪問者が式を訪れて記憶を得るメソッド `Q.visit`
* 変数の値のリストの下で、式訪問者が式を訪れて仕事をした結果を得る函数 `(::ExprVisitor)(::Expression, ::ValueList)`
* 代入する式訪問者の型 `Substitution`

```julia tags=[]
module R

########## 函数呼び出し機能を追加

using ..Q: Q, Expression

struct FunctionCall{F, X<:Expression} <: Expression f::F; x::X end
fun(fc::FunctionCall) = getfield(fc, :f)
arg(fc::FunctionCall) = getfield(fc, :x)

# 多重ディスパッチを利用
function Q.visit(visitor::Q.Stringifier, fc::FunctionCall)
    Q.visit(visitor, arg(fc))
    Q.memory(visitor)[fc] = string(fun(fc)) * "(" * Q.memory(visitor)[arg(fc)] * ")"
end

# 多重ディスパッチを利用
function Q.visit(visitor::Q.Evaluator, fc::FunctionCall)
    Q.visit(visitor, arg(fc))
    Q.memory(visitor)[fc] = fun(fc)(Q.memory(visitor)[arg(fc)])
end

########## 変数を扱う機能を追加

using ..Q: Q, Expression, ExprVisitor, Constant, BinOp, lhs, rhs

struct Variable <:Expression name::Symbol end
name(var::Variable) = getfield(var, :name)

# 多重ディスパッチを利用
function Q.visit(visitor::Q.Stringifier, var::Variable)
    Q.memory(visitor)[var] = string(name(var))
end

struct ValueList{T<:NamedTuple} v::T end
ValueList() = ValueList((;))
ValueList(; p...) = ValueList((; p...))
parent(list::ValueList) = getfield(list, :v)
names(list::ValueList) = keys(parent(list))
value(list::ValueList, name::Symbol) = getproperty(parent(list), name)
stringify(list::ValueList) = "ValueList" * string(parent(list))

# 多重ディスパッチを利用
function (visitor::ExprVisitor)(expr::Expression, list::ValueList)
    Q.visit(visitor, expr, list)
    Q.memory(visitor)[expr]
end

# 多重ディスパッチを利用
function Q.visit(visitor::Q.Evaluator, var::Variable, list::ValueList)
    Q.memory(visitor)[var] = value(list, name(var))
end
function Q.visit(visitor::Q.Evaluator, fc::FunctionCall, list::ValueList)
    Q.visit(visitor, arg(fc), list)
    Q.memory(visitor)[fc] = fun(fc)(Q.memory(visitor)[arg(fc)])
end
function Q.visit(visitor::Q.Evaluator, c::Constant, v::ValueList)
    Q.memory(visitor)[c] = Q.value(c)
end
for (S, op) in Q.binop_list
    @eval function Q.visit(visitor::Q.Evaluator, bo::Q.$S, list::ValueList=ValueList())
        Q.visit(visitor, lhs(bo), list)
        Q.visit(visitor, rhs(bo), list)
        Q.memory(visitor)[bo] = $op(Q.memory(visitor)[lhs(bo)], Q.memory(visitor)[rhs(bo)])
    end
end

struct Substitution <: ExprVisitor memory::Dict{Expression, Expression} end
Substitution() = Substitution(Dict{Expression, Expression}())

# 多重ディスパッチを利用
function Q.visit(visitor::Substitution, var::Variable, list::ValueList)
    Q.memory(visitor)[var] = name(var) ∈ names(list) ? Constant(value(list, name(var))) : var
end
function Q.visit(visitor::Substitution, fc::FunctionCall, list::ValueList)
    Q.visit(visitor, arg(fc), list)
    Q.memory(visitor)[fc] = FunctionCall(fun(fc), Q.memory(visitor)[arg(fc)])
end
function Q.visit(visitor::Substitution, c::Constant, list::ValueList)
    Q.memory(visitor)[c] = c
end
for (S, op) in Q.binop_list
    @eval function Q.visit(visitor::Substitution, bo::Q.$S, list::ValueList)
        Q.visit(visitor, lhs(bo), list)
        Q.visit(visitor, rhs(bo), list)
        Q.memory(visitor)[bo] = Q.$S(Q.memory(visitor)[lhs(bo)], Q.memory(visitor)[rhs(bo)])
    end
end

# デフォルトでの表示の仕方
Base.show(io::IO, list::ValueList) = print(io, stringify(list))

end
```

```julia
fc = R.FunctionCall(sinpi, Q.Div(Q.Constant(1), Q.Constant(6)))
```

```julia
stringifier = Q.Stringifier()
Q.visit(stringifier, fc)
Q.memory(stringifier)
```

```julia tags=[]
Q.Stringifier()(fc)
```

```julia
evaluator = Q.Evaluator()
Q.visit(evaluator, fc)
Q.memory(evaluator)
```

```julia tags=[]
Q.Evaluator()(fc)
```

```julia
fxy = R.FunctionCall(sinpi, Q.Div(R.Variable(:x), R.Variable(:y)))
```

```julia
substitution = R.Substitution()
Q.visit(substitution, fxy, R.ValueList(x = 1))
Q.memory(substitution)
```

```julia
R.Substitution()(fxy, R.ValueList(x = 1))
```

```julia
substitution = R.Substitution()
Q.visit(substitution, fxy, R.ValueList(x = 1, y = 6))
Q.memory(substitution)
```

```julia
R.Substitution()(fxy, R.ValueList(x = 1, y = 6))
```

```julia
Q.Evaluator()(fxy, R.ValueList(x = 1, y = 6))
```

```julia

```
