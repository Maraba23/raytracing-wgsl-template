const THREAD_COUNT = 16;
const RAY_TMIN = 0.0001;
const RAY_TMAX = 100.0;
const PI = 3.1415927f;
const FRAC_1_PI = 0.31830987f;
const FRAC_2_PI = 1.5707964f;

@group(0) @binding(0)  
  var<storage, read_write> fb : array<vec4f>;

@group(0) @binding(1)
  var<storage, read_write> rtfb : array<vec4f>;

@group(1) @binding(0)
  var<storage, read_write> uniforms : array<f32>;

@group(2) @binding(0)
  var<storage, read_write> spheresb : array<sphere>;

@group(2) @binding(1)
  var<storage, read_write> quadsb : array<quad>;

@group(2) @binding(2)
  var<storage, read_write> boxesb : array<box>;

@group(2) @binding(3)
  var<storage, read_write> trianglesb : array<triangle>;

@group(2) @binding(4)
  var<storage, read_write> meshb : array<mesh>;

struct ray {
  origin : vec3f,
  direction : vec3f,
};

struct sphere {
  transform : vec4f,
  color : vec4f,
  material : vec4f,
};

struct quad {
  Q : vec4f,
  u : vec4f,
  v : vec4f,
  color : vec4f,
  material : vec4f,
};

struct box {
  center : vec4f,
  radius : vec4f,
  rotation: vec4f,
  color : vec4f,
  material : vec4f,
};

struct triangle {
  v0 : vec4f,
  v1 : vec4f,
  v2 : vec4f,
};

struct mesh {
  transform : vec4f,
  scale : vec4f,
  rotation : vec4f,
  color : vec4f,
  material : vec4f,
  min : vec4f,
  max : vec4f,
  show_bb : f32,
  start : f32,
  end : f32,
};

struct material_behaviour {
  scatter : bool,
  direction : vec3f,
};

struct camera {
  origin : vec3f,
  lower_left_corner : vec3f,
  horizontal : vec3f,
  vertical : vec3f,
  u : vec3f,
  v : vec3f,
  w : vec3f,
  lens_radius : f32,
};

struct hit_record {
  t : f32,
  p : vec3f,
  normal : vec3f,
  object_color : vec4f,
  object_material : vec4f,
  frontface : bool,
  hit_anything : bool,
};

fn ray_at(r: ray, t: f32) -> vec3f
{
  return r.origin + t * r.direction;
}

fn get_ray(cam: camera, uv: vec2f, rng_state: ptr<function, u32>) -> ray
{
  var rd = cam.lens_radius * rng_next_vec3_in_unit_disk(rng_state);
  var offset = cam.u * rd.x + cam.v * rd.y;
  return ray(cam.origin + offset, normalize(cam.lower_left_corner + uv.x * cam.horizontal + uv.y * cam.vertical - cam.origin - offset));
}

fn get_camera(lookfrom: vec3f, lookat: vec3f, vup: vec3f, vfov: f32, aspect_ratio: f32, aperture: f32, focus_dist: f32) -> camera
{
  var camera = camera();
  camera.lens_radius = aperture / 2.0;

  var theta = degrees_to_radians(vfov);
  var h = tan(theta / 2.0);
  var w = aspect_ratio * h;

  camera.origin = lookfrom;
  camera.w = normalize(lookfrom - lookat);
  camera.u = normalize(cross(vup, camera.w));
  camera.v = cross(camera.u, camera.w);

  camera.lower_left_corner = camera.origin - w * focus_dist * camera.u - h * focus_dist * camera.v - focus_dist * camera.w;
  camera.horizontal = 2.0 * w * focus_dist * camera.u;
  camera.vertical = 2.0 * h * focus_dist * camera.v;

  return camera;
}

fn envoriment_color(direction: vec3f, color1: vec3f, color2: vec3f) -> vec3f
{
  var unit_direction = normalize(direction);
  var t = 0.5 * (unit_direction.y + 1.0);
  var col = (1.0 - t) * color1 + t * color2;

  var sun_direction = normalize(vec3(uniforms[13], uniforms[14], uniforms[15]));
  var sun_color = int_to_rgb(i32(uniforms[17]));
  var sun_intensity = uniforms[16];
  var sun_size = uniforms[18];

  var sun = clamp(dot(sun_direction, unit_direction), 0.0, 1.0);
  col += sun_color * max(0, (pow(sun, sun_size) * sun_intensity));

  return col;
}

// Função auxiliar para atualizar o registro mais próximo
fn update_closest(new_record: hit_record, object_color: vec4f, object_material: vec4f, closest: ptr<function, hit_record>) {
    if (new_record.hit_anything && new_record.t < (*closest).t) {
        *closest = new_record;
        (*closest).object_color = object_color;
        (*closest).object_material = object_material;
    }
}

fn check_ray_collision(r: ray, max: f32) -> hit_record {
    var spheresCount = i32(uniforms[19]);
    var quadsCount = i32(uniforms[20]);
    var boxesCount = i32(uniforms[21]);
    var meshCount = i32(uniforms[27]);

    var closest = hit_record(
        RAY_TMAX,
        vec3f(0.0),
        vec3f(0.0),
        vec4f(0.0),
        vec4f(0.0),
        false,
        false
    );
    var new_record = hit_record(
        RAY_TMAX,
        vec3f(0.0),
        vec3f(0.0),
        vec4f(0.0),
        vec4f(0.0),
        false,
        false
    );

    // Processa as esferas
    for (var i = 0; i < spheresCount; i++) {
        var s = spheresb[i];
        new_record = hit_record(
            RAY_TMAX,
            vec3f(0.0),
            vec3f(0.0),
            vec4f(0.0),
            vec4f(0.0),
            false,
            false
        );

        hit_sphere(s.transform.xyz, s.transform.w, r, &new_record, max);
        update_closest(new_record, s.color, s.material, &closest);
    }

    // Processa os quadrados
    for (var i = 0; i < quadsCount; i++) {
        var q = quadsb[i];
        new_record = hit_record(
            RAY_TMAX,
            vec3f(0.0),
            vec3f(0.0),
            vec4f(0.0),
            vec4f(0.0),
            false,
            false
        );

        hit_quad(r, q.Q, q.u, q.v, &new_record, max);
        update_closest(new_record, q.color, q.material, &closest);
    }

    // Processa as caixas
    for (var i = 0; i < boxesCount; i++) {
        var b = boxesb[i];
        new_record = hit_record(
            RAY_TMAX,
            vec3f(0.0),
            vec3f(0.0),
            vec4f(0.0),
            vec4f(0.0),
            false,
            false
        );

        hit_box(r, b.center.xyz, b.radius.xyz, &new_record, max);
        update_closest(new_record, b.color, b.material, &closest);
    }

    // Processa os meshes
    for (var i = 0; i < meshCount; i++) {
        var m = meshb[i];

        if (m.show_bb > 0.0) {
            new_record = hit_record(
                RAY_TMAX,
                vec3f(0.0),
                vec3f(0.0),
                vec4f(0.0),
                vec4f(0.0),
                false,
                false
            );

            var center = (m.min.xyz + m.max.xyz) * 0.5;
            var radius = (m.max.xyz - m.min.xyz) * 0.5;

            hit_box(r, center, radius, &new_record, max);
            update_closest(new_record, m.color, m.material, &closest);
        } else {
            for (var j = i32(m.start); j < i32(m.end); j++) {
                var t = trianglesb[j];
                new_record = hit_record(
                    RAY_TMAX,
                    vec3f(0.0),
                    vec3f(0.0),
                    vec4f(0.0),
                    vec4f(0.0),
                    false,
                    false
                );

                hit_triangle(r, t.v0.xyz, t.v1.xyz, t.v2.xyz, &new_record, max);
                update_closest(new_record, m.color, m.material, &closest);
            }
        }
    }

    return closest;
}


fn lambertian(normal : vec3f, absorption: f32, random_sphere: vec3f, rng_state: ptr<function, u32>) -> material_behaviour
{
  var res = normal + random_sphere;

  return material_behaviour(true, normalize(res));
}

fn metal(normal : vec3f, direction: vec3f, fuzz: f32, random_sphere: vec3f) -> material_behaviour
{
  var reflected = reflect(direction, normal);
  return material_behaviour(true, normalize(reflected + fuzz * random_sphere));
}

fn dielectric(
    normal: vec3f,
    r_direction: vec3f,
    refraction_index: f32,
    frontface: bool,
    random_sphere: vec3f,
    fuzz: f32,
    rng_state: ptr<function, u32>
) -> material_behaviour {
    // Determina a razão de índices de refração
    var eta_over_eta_prime: f32;
    if (frontface) {
        eta_over_eta_prime = 1.0 / refraction_index;
    } else {
        eta_over_eta_prime = refraction_index;
    }

    // Calcula cos_theta e sin_theta
    let cos_theta = dot(-r_direction, normal);
    let sin_theta = sqrt(1.0 - cos_theta * cos_theta);

    // Verifica se ocorre reflexão total interna
    let cannot_refract = eta_over_eta_prime * sin_theta > 1.0;

    // Calcula a refletância usando a aproximação de Schlick
    var reflectance = schlick_reflectance(cos_theta, refraction_index);

    // Gera um número aleatório para decidir entre refletir ou refratar
    let random_val = rng_next_float(rng_state);

    // Calcula a direção refletida
    let reflected = reflect(r_direction, normal);

    // Decide entre reflexão e refração
    var direction: vec3f;
    if (cannot_refract || random_val < reflectance) {
        direction = reflected;
    } else {
        // Calcula a direção refratada usando a Lei de Snell
        let r_out_perp = eta_over_eta_prime * (r_direction + cos_theta * normal);
        let r_out_parallel = -sqrt(abs(1.0 - dot(r_out_perp, r_out_perp))) * normal;
        direction = r_out_perp + r_out_parallel;
    }

    return material_behaviour(true, direction);
}

// Função auxiliar para calcular a refletância usando a aproximação de Schlick
fn schlick_reflectance(cos_theta: f32, refraction_index: f32) -> f32 {
    let r0 = pow((1.0 - refraction_index) / (1.0 + refraction_index), 2.0);
    return r0 + (1.0 - r0) * pow(1.0 - cos_theta, 5.0);
}


fn emmisive(color: vec3f, light: f32) -> material_behaviour
{
  return material_behaviour(false, light * color);
}

fn trace(r: ray, rng_state: ptr<function, u32>) -> vec3f
{
  var maxbounces = i32(uniforms[2]);
  var light = vec3f(0.0);
  var color = vec3f(1.0);
  var r_ = r;
  
  var backgroundcolor1 = int_to_rgb(i32(uniforms[11]));
  var backgroundcolor2 = int_to_rgb(i32(uniforms[12]));
  var behaviour = material_behaviour(true, vec3f(0.0));

  for (var j = 0; j < maxbounces; j = j + 1)
  {

    if (!behaviour.scatter) {
      break;
    }

    var closest_record = check_ray_collision(r_, RAY_TMAX);

    if (closest_record.hit_anything) {
      var normal = closest_record.normal;
      var prob = closest_record.object_material.z;
      var rand = rng_next_float(rng_state);

      var smoothness = closest_record.object_material.x * f32(rand < prob);

      var lamb_behaviour = lambertian(normal, closest_record.object_material.y, rng_next_vec3_in_unit_sphere(rng_state), rng_state); // WHAT ARE THE ARGS?
      var metal_behaviour = metal(normal, r_.direction, closest_record.object_material.y, rng_next_vec3_in_unit_sphere(rng_state)); // WHAT ARE THE ARGS?
      var emissive_behaviour = emmisive(closest_record.object_color.xyz, closest_record.object_material.w);
      var dielectric_behaviour = dielectric(normal, r_.direction, 1.0 + closest_record.object_material.z, closest_record.frontface, rng_next_vec3_in_unit_sphere(rng_state), closest_record.object_material.y, rng_state);

      var is_emissive = closest_record.object_material.w > 0;
      var is_dielectric = closest_record.object_material.x < 0;

      behaviour.scatter = lamb_behaviour.scatter && !is_emissive;

      if (is_dielectric) {
          behaviour.direction = dielectric_behaviour.direction;
      } else {
          behaviour.direction = mix(lamb_behaviour.direction, metal_behaviour.direction, smoothness);
      }

      var new_dir = behaviour.direction;
      var new_color = mix(closest_record.object_color.xyz, vec3(1.0), smoothness) * f32(!is_emissive) + emissive_behaviour.direction * f32(is_emissive);
      
      if (is_dielectric) {
          new_color = closest_record.object_color.xyz;
      }

      color *= new_color;

      r_ = ray(closest_record.p, new_dir);
      
      light += color * f32(is_emissive);
      
    } else {
      light += envoriment_color(r_.direction, backgroundcolor1, backgroundcolor2) * color;
      break;
    }
  }


  return light;
}

@compute @workgroup_size(THREAD_COUNT, THREAD_COUNT, 1)
fn render(@builtin(global_invocation_id) id: vec3u) {
    var resolution = uniforms[1];
    var frame_time = u32(uniforms[0]);
    var rng_state = init_rng(vec2(id.x, id.y), vec2(u32(resolution)), frame_time);

    var frag_coord = vec2f(f32(id.x), f32(id.y));
    var uv = (frag_coord + sample_square(&rng_state)) / vec2(resolution);

    var camera_position = vec3(uniforms[7], uniforms[8], uniforms[9]);
    var target_position = vec3(uniforms[23], uniforms[24], uniforms[25]);
    var camera = get_camera(
        camera_position, 
        target_position, 
        vec3(0.0, 1.0, 0.0), 
        uniforms[10], 
        1.0, 
        uniforms[6], 
        uniforms[5]
    );

    var samples_per_pixel = i32(uniforms[4]);
    var pixel_color = vec3(0.0);

    for (var sample = 0; sample < samples_per_pixel; sample++) {
        var ray = get_ray(camera, uv, &rng_state);
        pixel_color += trace(ray, &rng_state);
    }

    pixel_color /= f32(samples_per_pixel);
    pixel_color = saturate(pixel_color);

    var final_color = vec4(linear_to_gamma(pixel_color), 1.0);
    var framebuffer_index = mapfb(id.xy, resolution);
    var should_accumulate = uniforms[3];
    var accumulated_color = rtfb[framebuffer_index] * should_accumulate + final_color;
    var weight = accumulated_color.w;

    rtfb[framebuffer_index] = accumulated_color;
    fb[framebuffer_index] = accumulated_color / weight;
}
