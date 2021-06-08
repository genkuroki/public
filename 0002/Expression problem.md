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
value(x::Constant) = getfield(x, :value)
stringify(x::Constant) = string(value(x))
evaluate(x::Constant) = value(x)

abstract type BinOp <: Expression end
lhs(x::BinOp) = getfield(x, :lhs)
rhs(x::BinOp) = getfield(x, :rhs)

const binop_list = ((:Plus, :+), (:Minus, :-), (:Mul, :*), (:Div, :/))
for (S, op) in binop_list
    @eval begin
        struct $S{L<:Expression, R<:Expression} <: BinOp lhs::L; rhs::R end
        stringify(x::$S) = "(" * stringify(lhs(x)) * " $($op) " * stringify(rhs(x)) * ")"
        evaluate(x::$S) = $op(evaluate(lhs(x)), evaluate(rhs(x)))
    end
end

# デフォルトでの表示の仕方
Base.show(io::IO, x::Expression) = print(io, stringify(x))

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
fun(fx::FunctionCall) = getfield(fx, :f)
arg(fx::FunctionCall) = getfield(fx, :x)
O.stringify(fx::FunctionCall) = string(fun(fx)) * "(" * O.stringify(arg(fx)) * ")"
O.evaluate(fx::FunctionCall) = fun(fx)(O.evaluate(arg(fx)))

### 変数を扱う機能を追加

using ..O: O, Expression, Constant, BinOp, lhs, rhs
using InteractiveUtils: subtypes

struct Variable <:Expression name::Symbol end
name(x::Variable) = getfield(x, :name)
O.stringify(x::Variable) = string(name(x))

struct ValueList{T<:NamedTuple} v::T end
ValueList() = ValueList((;))
ValueList(; p...) = ValueList((; p...))
parent(x::ValueList) = getfield(x, :v)
names(x::ValueList) = keys(parent(x))
value(x::ValueList, name::Symbol) = getproperty(parent(x), name)
O.stringify(x::ValueList) = "ValueList" * string(parent(x))

evaluate(x::Variable, v::ValueList) = value(v, name(x))
evaluate(x::FunctionCall, v::ValueList) = fun(x)(evaluate(arg(x), v))
evaluate(x::Constant, v::ValueList) = O.value(x)
for (S, op) in O.binop_list
    @eval evaluate(x::O.$S, v::ValueList=ValueList()) = $op(evaluate(lhs(x), v), evaluate(rhs(x), v))
end

substitute(x::Variable, v::ValueList) = name(x) ∈ names(v) ? Constant(value(v, name(x))) : x
substitute(x::FunctionCall, v::ValueList) = FunctionCall(fun(x), substitute(arg(x), v))
substitute(x::Constant, v::ValueList) = x
for (S, op) in O.binop_list
    @eval substitute(x::O.$S, v::ValueList) = O.$S(substitute(lhs(x), v), substitute(rhs(x), v))
end

# デフォルトでの表示の仕方
Base.show(io::IO, x::ValueList) = print(io, O.stringify(x))

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

* https://eli.thegreenplace.net/2016/the-expression-problem-and-its-solutions/

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

Visitor patternに従って、モジュールOを完全に作り直してしまう。 新モジュール名をQとする。

モジュールOと比較すると圧倒的に複雑で面倒になってしまっている。

* 式の型 `Expression`
  * 定数の型 `Constant`
    * 定数の値を得る函数 `value`
  * 二項演算の型 `Plus`, `Minus`, `Mul`, `Div`
    * 二項演算の左辺と右辺を得る函数 `lhs`, `rhs`
* 式訪問者の型 `ExprVisitor`
  * 式を文字列化する者の型 `Stringifier`
  * 式の評価者の型 `Evaluator`
  * 式訪問者の記憶を取り出す函数 `memory`
  * 式訪問者が式を訪れて仕事をした結果を得る函数 `(::ExprVisitor)(::Expression)`
  * 式訪問者が式を訪れて記憶(思い出)を得る函数 `visit`

```julia
module Q

abstract type Expression end

struct Constant{T} <: Expression value::T end
value(x::Constant) = getfield(x, :value)

abstract type BinOp <: Expression end
lhs(x::BinOp) = getfield(x, :lhs)
rhs(x::BinOp) = getfield(x, :rhs)

const binop_list = ((:Plus, :+), (:Minus, :-), (:Mul, :*), (:Div, :/))
for (S, op) in binop_list
    @eval struct $S{L<:Expression, R<:Expression} <: BinOp lhs::L; rhs::R end
end

abstract type ExprVisitor end
memory(x::ExprVisitor) = getfield(x, :memory)
function (x::ExprVisitor)(e::Expression) # 多重ディスパッチを利用
    visit(x, e)
    memory(x)[e]
end

struct Evaluator <: ExprVisitor memory::Dict{Expression, Any} end
Evaluator() = Evaluator(Dict{Expression, Any}())

# 多重ディスパッチを利用
visit(x::Evaluator, c::Constant) = memory(x)[c] = value(c)
for (S, op) in binop_list
    @eval function visit(x::Evaluator, bp::$S)
        visit(x, lhs(bp))
        visit(x, rhs(bp))
        memory(x)[bp] = $op(memory(x)[lhs(bp)], memory(x)[rhs(bp)])
    end
end

struct Stringifier <: ExprVisitor memory::Dict{Expression, String} end
Stringifier() = Stringifier(Dict{Expression, String}())

# 多重ディスパッチを利用
visit(x::Stringifier, c::Constant) = memory(x)[c] = string(value(c))
for (S, op) in binop_list
    @eval function visit(x::Stringifier, bp::$S)
        visit(x, lhs(bp))
        visit(x, rhs(bp))
        memory(x)[bp] = "(" * memory(x)[lhs(bp)] * " $($op) " * memory(x)[rhs(bp)] * ")"
    end
end

# デフォルトでの表示の仕方
Base.show(io::IO, x::Expression) = print(io, Stringifier()(x))

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
## Visitor patternで機能を追加

* https://eli.thegreenplace.net/2016/the-expression-problem-and-its-solutions/

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

多重ディスパッチのおかげで、上のＣ＋＋版と比較すると圧倒的にシンプルになっている。

しかし、モジュールPにおける次のコードと比較するとロジックは複雑で込み入ったものになってしまう。

```Julia
struct FunctionCall{F, X<:Expression} <: Expression f::F; x::X end
fun(fx::FunctionCall) = getfield(fx, :f)
arg(fx::FunctionCall) = getfield(fx, :x)
O.stringify(fx::FunctionCall) = string(fun(fx)) * "(" * O.stringify(arg(fx)) * ")"
O.evaluate(fx::FunctionCall) = fun(fx)(O.evaluate(arg(fx)))
```

以下を追加する。

* 函数呼び出しの型 `FunctionCall`
  * 函数呼び出しの函数部分を取り出す函数 `fun`
  * 函数呼び出しの引数部分を取り出す函数 `arg`
* 式を文字列に変換する者が函数呼び出しを訪れて記憶を得るメソッド `Q.visit`
* 式評価者が函数呼び出しを訪れて記憶を得るメソッド `Q.visit`

```julia tags=[]
module R

### 函数呼び出し機能を追加

using ..Q: Q, Expression

struct FunctionCall{F, X<:Expression} <: Expression f::F; x::X end
fun(fx::FunctionCall) = getfield(fx, :f)
arg(fx::FunctionCall) = getfield(fx, :x)

# 多重ディスパッチを利用
function Q.visit(x::Q.Stringifier, fx::FunctionCall)
    Q.visit(x, arg(fx))
    Q.memory(x)[fx] = string(fun(fx)) * "(" * Q.memory(x)[arg(fx)] * ")"
end

# 多重ディスパッチを利用
function Q.visit(x::Q.Evaluator, fx::FunctionCall)
    Q.visit(x, arg(fx))
    Q.memory(x)[fx] = fun(fx)(Q.memory(x)[arg(fx)])
end

end
```

```julia
fx = R.FunctionCall(sinpi, Q.Div(Q.Constant(1), Q.Constant(6)))
```

```julia
stringifier = Q.Stringifier()
Q.visit(stringifier, fx)
Q.memory(stringifier)
```

```julia tags=[]
Q.Stringifier()(fx)
```

```julia
evaluator = Q.Evaluator()
Q.visit(evaluator, fx)
Q.memory(evaluator)
```

```julia tags=[]
Q.Evaluator()(fx)
```

```julia

```
