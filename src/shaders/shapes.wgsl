fn hit_sphere(center: vec3f, radius: f32, r: ray, record: ptr<function, hit_record>, max: f32) {
    // Vetor do centro da esfera até a origem do raio
    let oc = r.origin - center;
    let a = dot(r.direction, r.direction);
    let half_b = dot(oc, r.direction);
    let c = dot(oc, oc) - radius * radius;
    let discriminant = half_b * half_b - a * c;

    if (discriminant < 0.0) {
        (*record).hit_anything = false;
        return;
    }

    let sqrt_d = sqrt(discriminant);

    // Encontra a raiz mais próxima que está dentro do intervalo válido
    var t = (-half_b - sqrt_d) / a;
    if (t < RAY_TMIN || t > max) {
        t = (-half_b + sqrt_d) / a;
        if (t < RAY_TMIN || t > max) {
            (*record).hit_anything = false;
            return;
        }
    }

    // Atualiza o registro de interseção
    let p = ray_at(r, t);
    let outward_normal = normalize(p - center);
    (*record).t = t;
    (*record).p = p;

    // Determina se o raio está atingindo a frente ou o verso da superfície
    if (dot(r.direction, outward_normal) < 0.0) {
        (*record).frontface = true;
        (*record).normal = outward_normal;
    } else {
        (*record).frontface = false;
        (*record).normal = -outward_normal;
    }

    (*record).hit_anything = true;
}


fn hit_quad(r: ray, Q: vec4f, u: vec4f, v: vec4f, record: ptr<function, hit_record>, max: f32) {
    // Calcula o vetor normal do plano do quadrado
    let n = cross(u.xyz, v.xyz);
    let normal = normalize(n);

    // Calcula o denominador para a equação de interseção
    let denom = dot(normal, r.direction);

    // Verifica se o raio é paralelo ao plano
    if abs(denom) < 1e-4 {
        (*record).hit_anything = false;
        return;
    }

    // Calcula o parâmetro t para a interseção
    let t = dot(normal, Q.xyz - r.origin) / denom;

    // Verifica se t está dentro do intervalo válido
    if t < RAY_TMIN || t > max {
        (*record).hit_anything = false;
        return;
    }

    // Calcula o ponto de interseção
    let p = ray_at(r, t);

    // Vetor do vértice Q até o ponto de interseção
    let d = p - Q.xyz;

    // Calcula os produtos escalares necessários
    let ddotu = dot(d, u.xyz);
    let ddotv = dot(d, v.xyz);
    let udotu = dot(u.xyz, u.xyz);
    let vdotv = dot(v.xyz, v.xyz);

    // Calcula os parâmetros alfa e beta
    let alpha = ddotu / udotu;
    let beta = ddotv / vdotv;

    // Verifica se o ponto está dentro do quadrado
    if alpha < 0.0 || alpha > 1.0 || beta < 0.0 || beta > 1.0 {
        (*record).hit_anything = false;
        return;
    }

    // Determina se o raio está atingindo a frente ou o verso do quadrado
    let frontface = dot(normal, r.direction) < 0.0;

    // Ignora interseções com a face posterior do quadrado
    if !frontface {
        (*record).hit_anything = false;
        return;
    }

    // Atualiza o registro de interseção
    (*record).t = t;
    (*record).p = p;
    (*record).frontface = frontface;
    (*record).normal = normal;
    (*record).hit_anything = true;
}


fn hit_triangle(r: ray, v0: vec3f, v1: vec3f, v2: vec3f, record: ptr<function, hit_record>, max: f32) {
    // Vetores das arestas do triângulo
    let edge1 = v1 - v0;
    let edge2 = v2 - v0;

    // Vetor normal do triângulo
    let h = cross(r.direction, edge2);
    let a = dot(edge1, h);

    // Verifica se o raio é paralelo ao triângulo
    if abs(a) < 1e-4 {
        (*record).hit_anything = false;
        return;
    }

    let f = 1.0 / a;
    let s = r.origin - v0;
    let u = f * dot(s, h);

    // Verifica se o ponto está fora do triângulo
    if u < 0.0 || u > 1.0 {
        (*record).hit_anything = false;
        return;
    }

    let q = cross(s, edge1);
    let v = f * dot(r.direction, q);

    if v < 0.0 || (u + v) > 1.0 {
        (*record).hit_anything = false;
        return;
    }

    // Calcula t para encontrar o ponto de interseção ao longo do raio
    let t = f * dot(edge2, q);

    if t < RAY_TMIN || t > max {
        (*record).hit_anything = false;
        return;
    }

    // Atualiza o registro de interseção
    (*record).t = t;
    (*record).p = ray_at(r, t);
    let normal = normalize(cross(edge1, edge2));

    // Determina se o raio está atingindo a frente ou o verso do triângulo
    let frontface = dot(r.direction, normal) < 0.0;

    // Ajusta o normal com base na face atingida
    if frontface {
        (*record).normal = normal;
    } else {
        (*record).normal = -normal;
    }

    (*record).frontface = frontface;
    (*record).hit_anything = true;
}


fn hit_box(r: ray, center: vec3f, rad: vec3f, record: ptr<function, hit_record>, t_max: f32)
{
  var m = 1.0 / r.direction;
  var n = m * (r.origin - center);
  var k = abs(m) * rad;

  var t1 = -n - k;
  var t2 = -n + k;

  var tN = max(max(t1.x, t1.y), t1.z);
  var tF = min(min(t2.x, t2.y), t2.z);

  if (tN > tF || tF < 0.0)
  {
    (*record).hit_anything = false;
    return;
  }

  var t = tN;
  if (t < RAY_TMIN || t > t_max)
  {
    (*record).hit_anything = false;
    return;
  }

  (*record).t = t;
  (*record).p = ray_at(r, t);
  (*record).normal = -sign(r.direction) * step(t1.yzx, t1.xyz) * step(t1.zxy, t1.xyz);
  (*record).hit_anything = true;

  return;
}
