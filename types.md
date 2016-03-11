### Some things with mypy: sum and product types
[A recent pep](https://www.python.org/dev/peps/pep-0484/) solidifies type annotations in python 2 and 3. These type annotations are compatible with current versions of python. 
They can be ignored altogether, used however you wish in your own program, or used to typecheck your code with [mypy](github.com/python/mypyp).
This post will discuss the last option. Later we'll see how python's strong introspective powers can be used to leverage these types in conjunction with
mypy.

`mypy` enables static typing in python. The features include defacto case-classes (using named tuples) and Union types. These are known
more generally as "product types" and "sum types" respectively. A product type is similar to a class in Java. It has pre-defined members (private or public)
of other types. In a sense it is a "product" of these other types. In mypy, one can declare a product type using classes, or more simply using `NamedTuple.`
For example, let's create a product type for points in a 3D plane.
```python
Point3D = NamedTuple("Point3D", [("x", float), ("y", float), ("z", float)])
```
If we wanted to use a simple tuple instead, we could declare that so:

```python
Point3DTuple = Tuple[float, float, float]
```

Let's look at what our named tuple can, and more importantly can't, do within mypy's type world. All the displayed errors
are part of mypy's output, which provides the line of the error as well as an explanation. Keep in mind these are all erorrs caught before
the program, or even any tests, are run. These errors can also be integrated with editors or IDEs to provde real-time feedback.

It can be created and accessed just like `collections.namedtuple.`
```python
point = Point3D(0, 1.0, 3.98)
x = point.x
y = point[1] # this typechecks, it probably shouldn't
```
mypy knows how long the tuple is, and what types its members are!
```python
r = point[99]
foo.py:10: error: Tuple index out of range
```
mypy enforces the safety of common operators. This avoids meaningless comparisons, for example, which are uncaught by python's runtime:
```python
>>> "foo" > sys.maxint
 True # sure, why not?
```
```python
point.x + "Eureka"
foo.py:10: error: Unsupported operand types for + ("float" and "str")

x = point.x # mypy infers the type after assignment

x > "Eureka"
foo.py:10: error: Unsupported operand types for > ("float" and "str")
```
mypy limits attribute access:

```python
sneaky = point.gecko
foo.py:13: error: "Point3D" has no attribute "gecko"
```
mypy supports generics. A generic can be a lot of things; A `list`, an `Iterable`, or something equivalent to scala/java 8's `Option` type. If a generic is a collection, all elements of the collection must be of the same type. mypy comes equipped with a number of generic types; take for example `List`, which is an alias for the built-in `list`.
```python 
ListOfInts = List[int]
```

You can also create types by subclassing `Generic`.
```python
class Option(Generic[T]):
    def getOrElse(t: T) -> T:
       . . . 
```
It's possible to use multiple type variables within a generic:
```python
E = TypeVar("E")
V = TypeVar("V")
class Either(Generic[E,V]):
    . . . . 
```

Let's use `List` and `3DPoint` to create a more complex product type: `Robot Legs`.

```python
RobotLegs = NamedTuple("RobotLegs", [("leftLeg", List[Point3D]), ("rightLeg", List[Point3D]), ("color", str)])
```
Note that we've defined the field `color` as simply a string, allowing us to create robot legs with nonsense colors. It's also possible to create robot legs with negative integers for coordinates! We only want pastel colors, and robots which exist in the cartesian plane. 
```python
blueRobot = RobotLegs(points, points, "fizbizzle")
```
Of course, we could check for this condition in the functions that use the color:
```python
def getColor(legs: RobotLegs) -> int:
    if legs.color not in ["skyblue", "red", "white"]:
        raise ValueError("Invalid color %s" % legs.color)
    else:
         . . . . 
```
That's a hassle, and it's easy to forget to do these checks in every function. Instead, let's nip this in the bud. 
We really want to make it is easy on ourselves and be *really really sure* that we only have to validate our input once. We can do all the validation--cleaning up data from I/O, verifying it matches a certain shape, creating errors etc.--when we construct the instances of our types. That way all functions which accept those types are relieved from the obligation of checking themselves.
```python
SkyBlue = NamedTuple("SkyBlue", [])
PastelRed = NamedTuple("PastelRed", [])
White = NamedTuple("White", [])

Color = Union[Blue, PastelRed, White]

RobotLegs = NamedTuple("RobotLegs", [("leftLeg", List[Point3D]), ("rightLeg", List[Point3D]), ("color", Color)])

make3DCoordinates(x: float, y: float, z: float) -> Point3D:
    assert x >= 0 and y >= 0 and z >= 0
    return Point3D(x, y, z)
```
Now we can be assured that our color is one of the primaries (always a good starting pint for giant robots), so we don't have to worry about validating our data again!

```python
def getColor(legs: RobotLegs) -> int:
    if legs.color == SkyBlue():  return 0x87CEFA 
    if isinstance(legs.color, SkyBlue): return  0x87CEFA # this is equivalent
```

We can even safely use a statically typed dictionary which never raise a KeyErorr:
```python
colors = { SkyBlue() : 0x87CEFA } # type: Dict[Color,int]
. . . . 
```

In fact it's possible to use this technique to *guarantee* that our function will only ever get valid input. It's only possible to construct the sum type of `RobotLegs` through the union type of `Color`; `Color` is by definition one of `Blue`, `Red`. . . and points
In languages with the concept of private constructors, it's possible to *guarantee* that a RobotLegs cannot be created an invalid state--and therefore that `getColor` can never be passed invalid data--by making the `RobotLegs` constructor private. Unfortunately, we can only document the `make3DCoordinates` function as the point of entry for our API--we can't exclude the constructor as private.

Note that the assurance offered by static typing is significantly stronger than the contract offered by ducked typing. If we simply accepted an object with `leftLeg` `rightLeg` and `color` as a RobotLeg, we'd have no guarantees that these fields were valid, or even that they were the expected type!

`Color` is a very simple Union type, analogous to the "Enums" of other languages (including python 3), while providing additional safety. Bution union types are more powerful; it's possible to create a union type out of product types, and model arbitrary complex 
systems this way. You can think of these types as representing the "set of all possible inputs and outputs" and functions accepting these types as representing the "cobminators" or "all the things I can ever do with my inputs". Together, these form a sort of "algebra" that represents your domain. In the domain of giant robots:

```python
Rifle = NamedTuple('Rifle', [('ammo' , int), ('model' , str)])
Knife = NamedTuple('Knife', [('shape' , List[Point3D]), ('thatsNotAKnife', bool)])

weapon = Union[Rifle, Knife]

RobotLegs = NamedTuple("RobotArms", [("leftArm", List[Point3D]), ("rightArm", List[Point3D]), ("color", Color)])

GiantRobot = NamedTuple('GiantRobot', [('weapon', Weapon), ('legs' , RobotLegs), ('arms', RobotArms)])

def canFight(robot: GiantRobot) -> bool:
    if isinstance(robot.weapon, Rifle):
        return robot.weapon.ammo > 0
    else: 
        return not robot.weapon.thatsNotAKnife # this is a knife
```
The `isinstance` check tells mypy that `robot.weapon` is specifically a rifle, and therefore has an `ammo` field of type `int`. Without that check, we get a nifty error from mypy--and find out before testing, running, or deploying:
```python
foo.py: note: In function "canFight":
foo.py:35: error: Some element of union has no attribute "ammo"
```
Great! we've created an API that's clear, self-documenting, and compartively safe. We've provided some limited guarantees of correctness;
and our domain is well-defined, which will help us reason about our past and future code moving forward.
mypy is a growing project; it's still in an early stage and being actively developed. It's become an official
part of they [python](github.com/python) flock as the definitive optional typechecker; it's got the [backing](https://github.com/python/mypy/issues/1276#issuecomment-192981427)
and [involvement](https://github.com/python/mypy/pull/1277) of [python's creator](https://en.wikipedia.org/wiki/Guido_van_Rossum).

Although mypy is still in active development, it can be a profitable tool right now. It's not a compiler, and it never touches
your code, so it can be used without much concern for bugs. It takes some extra time to annotate python with types--I've demonstrated
some of the strengths of its type inference, but it's necessary to annotate some things like lambda expressions, for example.
It's well worth the effort to document and verify your code in one way or another--mypy is another excellent tool for this purpose.
