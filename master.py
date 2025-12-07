import math
from typing import List, Tuple, Optional
import matplotlib.pyplot as plt
import matplotlib.patches as patches


class ConvexPolygon:
    def __init__(self, vertices: List[Tuple[float, float]]):

        if len(vertices) < 3:
            raise ValueError("многоугольник должен иметь хотя бы 3 вершины")

        self.vertices = vertices
        self.n = len(vertices)

        if not self._is_convex():
            raise ValueError("многоугольник не является выпуклым")

    def _cross_product(self, o: Tuple[float, float], a: Tuple[float, float], b: Tuple[float, float]) -> float:
        #векторное произведение (OAₓ × OBᵧ) - (OAᵧ × OBₓ)
        #если результат > 0 - поворот от OA к OB против часовой стрелки
        #если результат < 0 - поворот по часовой стрелке
        #если результат = 0 - точки коллинеарны
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

    def _is_convex(self) -> bool:
        if self.n < 3:
            return False

        # проверяем знак векторных произведений
        sign = None
        for i in range(self.n):
            o = self.vertices[i]
            a = self.vertices[(i + 1) % self.n]
            b = self.vertices[(i + 2) % self.n]

            cross = self._cross_product(o, a, b)

            # тут проверка на ноль, но пришлось поменять логику проверки на ноль из-за погрешности
            # плавающей точки. Число 1e-10 (0.0000000001) - это порог точности для учета ошибок округления
            # при работе с числами с плавающей точкой.
            if abs(cross) < 1e-10:  # коллинеарные точки
                continue

            # если у какой-то тройки знак отличается - многоугольник невыпуклый
            if sign is None:
                sign = cross > 0
            else:
                if (cross > 0) != sign:
                    return False

        return True

    # вычисление площади многоугольника методом шнуровки
    def area(self) -> float:
        # вычисляется путем суммирования произведений координат соседних вершин,
        # а затем нахождения половины модуля полученной суммы
        area = 0.0
        for i in range(self.n):
            j = (i + 1) % self.n
            area += self.vertices[i][0] * self.vertices[j][1]
            area -= self.vertices[j][0] * self.vertices[i][1]
        return abs(area) / 2.0

    # вычисление периметра многоугольника
    def perimeter(self) -> float:
        perimeter = 0.0
        for i in range(self.n):
            j = (i + 1) % self.n
            dx = self.vertices[j][0] - self.vertices[i][0]
            dy = self.vertices[j][1] - self.vertices[i][1]
            perimeter += math.sqrt(dx * dx + dy * dy)
        return perimeter

    # проверка, находится ли точка внутри многоугольника.
    def contains_point(self, point: Tuple[float, float]) -> bool:
        # используется метод лучей (проверка количества пересечений).
        x, y = point
        inside = False

        for i in range(self.n):
            j = (i + 1) % self.n
            xi, yi = self.vertices[i]
            xj, yj = self.vertices[j]

            if ((yi > y) != (yj > y)) and (x < (xj - xi) * (y - yi) / (yj - yi) + xi):
                inside = not inside

        return inside

    # проверка, находится ли другой многоугольник внутри данного
    def contains_polygon(self, other: 'ConvexPolygon') -> bool:
        return all(self.contains_point(vertex) for vertex in other.vertices)

    # нахождение точки пересечения двух отрезков
    def _line_intersection(self, p1: Tuple[float, float], p2: Tuple[float, float],
                           p3: Tuple[float, float], p4: Tuple[float, float]) -> Optional[Tuple[float, float]]:
        x1, y1 = p1
        x2, y2 = p2
        x3, y3 = p3
        x4, y4 = p4

        denom = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4)

        if abs(denom) < 1e-10:
            return None  # параллельные или совпадающие прямые

        t = ((x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4)) / denom
        u = -((x1 - x2) * (y1 - y3) - (y1 - y2) * (x1 - x3)) / denom

        if 0 <= t <= 1 and 0 <= u <= 1:
            x = x1 + t * (x2 - x1)
            y = y1 + t * (y2 - y1)
            return (x, y)

        return None


    # нахождение пересечения двух выпуклых многоугольников.
    def intersection(self, other: 'ConvexPolygon') -> Optional['ConvexPolygon']:
        # используется алгоритм Sutherland-Hodgman.
        # изначально использовалась проверка на принадлежность вершин, однако не учитываются
        # все случаи (звезда Давида)

        # Правки: мы избегаем лишних копирований, используя два списка
        # (current и next_poly) и переключаясь между ними.

        # проверка, находится ли точка p внутри ребра ab
        def inside(p: Tuple[float, float], a: Tuple[float, float], b: Tuple[float, float]) -> bool:
            return self._cross_product(a, b, p) >= 0

        # вычисление точки пересечения отрезка p1p2 с ребром ab
        def compute_intersection(p1: Tuple[float, float], p2: Tuple[float, float],
                                 a: Tuple[float, float], b: Tuple[float, float]) -> Tuple[float, float]:
            # уравнения прямых в общем виде: A*x + B*y = C
            A1 = p2[1] - p1[1]
            B1 = p1[0] - p2[0]
            C1 = A1 * p1[0] + B1 * p1[1]

            A2 = b[1] - a[1]
            B2 = a[0] - b[0]
            C2 = A2 * a[0] + B2 * a[1]

            det = A1 * B2 - A2 * B1  # определитель

            if abs(det) < 1e-10:
                return p1  # параллельные прямые

            x = (B2 * C1 - B1 * C2) / det
            y = (A1 * C2 - A2 * C1) / det

            return (x, y)

        # Используем два списка вместо копирования
        current = list(other.vertices)  # текущий многоугольник после обработки предыдущих ребер
        next_poly = []  # многоугольник после обработки текущим ребром

        # для каждого ребра clipping многоугольника
        for i in range(self.n):
            if not current:  # если многоугольник стал пустым
                break

            a = self.vertices[i]
            b = self.vertices[(i + 1) % self.n]

            # Очищаем список для нового многоугольника
            # Вместо создания нового списка, очищаем существующий
            next_poly.clear()

            # обходим все ребра текущего многоугольника
            m = len(current)
            for j in range(m):
                p1 = current[j]
                p2 = current[(j + 1) % m]

                # сохраняем только точки, лежащие внутри текущего ребра отсечения
                if inside(p2, a, b):
                    if not inside(p1, a, b):
                        next_poly.append(compute_intersection(p1, p2, a, b))
                    next_poly.append(p2)
                elif inside(p1, a, b):
                    next_poly.append(compute_intersection(p1, p2, a, b))

            # Переключаем списки: результат этой итерации становится входом для следующей
            current, next_poly = next_poly, current

        if len(current) < 3:
            return None

        try:
            return ConvexPolygon(current)
        except ValueError:
            # В случае численных ошибок (почти коллинеарные точки)
            return None

    # триангуляция выпуклого многоугольника.
    def triangulate(self) -> List[List[Tuple[float, float]]]:
        triangles = []

        # для выпуклого многоугольника простая триангуляция - веер из первой вершины
        # берем первую вершину как общую для всех треугольников
        for i in range(1, self.n - 1):
            triangle = [self.vertices[0], self.vertices[i], self.vertices[i + 1]]
            triangles.append(triangle) # для n-угольника получаем n-2 треугольника

        return triangles

    def plot(self, color='blue', alpha=0.5, label=None):
        polygon = patches.Polygon(self.vertices, closed=True,
                                  color=color, alpha=alpha, label=label)
        plt.gca().add_patch(polygon)

        x, y = zip(*self.vertices + [self.vertices[0]])
        plt.plot(x, y, 'o-', color=color, markersize=4)

    def __repr__(self):
        return f"ConvexPolygon({self.vertices})"


if __name__ == "__main__":
    square = ConvexPolygon([(0, 0), (2, 0), (2, 2), (0, 2)])

    triangle = ConvexPolygon([(1, 1), (3, 1), (2, 3)])

    pentagon = ConvexPolygon([(0.5, 0.5), (1.5, 0.2), (2.5, 0.5),
                              (2.2, 1.5), (0.8, 1.5)])

    print(f"площадь квадрата: {square.area():.2f}")
    print(f"периметр квадрата: {square.perimeter():.2f}")
    print(f"площадь треугольника: {triangle.area():.2f}")
    print(f"периметр треугольника: {triangle.perimeter():.2f}")

    test_point = (1, 1)
    print(f"точка {test_point} в квадрате: {square.contains_point(test_point)}")
    print(f"точка {test_point} в треугольнике: {triangle.contains_point(test_point)}")

    print(f"треугольник в квадрате: {square.contains_polygon(triangle)}")

    intersection = square.intersection(triangle)
    print(f"\nПересечение существует: {intersection is not None}")
    if intersection:
        print(f"Вершины пересечения: {intersection.vertices}")
        print(f"Площадь пересечения: {intersection.area():.2f}")
    else:
        print("Пересечение пустое или None")

    triangles = pentagon.triangulate()
    print(f"пятиугольник разбит на {len(triangles)} треугольника")

    plt.figure(figsize=(12, 10))

    plt.subplot(2, 2, 1)
    square.plot(color='blue', alpha=0.5, label='Квадрат')
    triangle.plot(color='red', alpha=0.5, label='Треугольник')
    plt.plot(test_point[0], test_point[1], 'go', markersize=8, label='Тестовая точка')
    plt.legend()
    plt.title('Исходные многоугольники')
    plt.axis('equal')
    plt.grid(True)

    plt.subplot(2, 2, 2)
    square.plot(color='blue', alpha=0.3, label='Квадрат')
    triangle.plot(color='red', alpha=0.3, label='Треугольник')
    if intersection:
        intersection.plot(color='green', alpha=0.7, label='Пересечение')
    plt.legend()
    plt.title('Пересечение многоугольников')
    plt.axis('equal')
    plt.grid(True)

    plt.subplot(2, 2, 3)
    colors = ['red', 'green', 'blue', 'orange']
    for i, triangle_vertices in enumerate(triangles):
        tri_patch = patches.Polygon(triangle_vertices, closed=True,
                                    color=colors[i % len(colors)], alpha=0.6)
        plt.gca().add_patch(tri_patch)

        x, y = zip(*triangle_vertices + [triangle_vertices[0]])
        plt.plot(x, y, 'o-', color=colors[i % len(colors)], markersize=4)

    plt.title('Триангуляция пятиугольника')
    plt.axis('equal')
    plt.grid(True)

    plt.subplot(2, 2, 4)

    convex_demo = ConvexPolygon([(0, 0), (2, 0), (2, 2), (0, 2)])

    convex_patch = patches.Polygon(convex_demo.vertices, closed=True,
                                   color='green', alpha=0.7, label='Выпуклый')
    plt.gca().add_patch(convex_patch)

    x1, y1 = zip(*convex_demo.vertices + [convex_demo.vertices[0]])
    plt.plot(x1, y1, 'g-', alpha=0.7)

    # пытаемся создать невыпуклый
    try:
        non_convex = ConvexPolygon([(0, 0), (2, 0), (1, 1), (2, 2), (0, 2)])
        non_convex.plot(color='red', alpha=0.7, label='Невыпуклый')
    except ValueError:
        non_convex_vertices = [(0, 0), (2, 0), (1, 1), (2, 2), (0, 2)]
        x2, y2 = zip(*non_convex_vertices + [non_convex_vertices[0]])
        plt.plot(x2, y2, 'ro--', alpha=0.7, linewidth=2, label='Невыпуклый (ошибка)')

    plt.legend()
    plt.title('Сравнение: выпуклый vs невыпуклый')
    plt.axis('equal')
    plt.grid(True)

    plt.tight_layout()
    plt.show()
